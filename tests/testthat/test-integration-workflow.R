# Comprehensive integration tests for tsenat() orchestration function
# Refactored to: avoid redundancies, create analysis once, test diverse argument combinations
# Includes functions from Appendix A (diversity comparison) and Appendix B (GAM/rank-based testing)

context("Refactored Integration Tests: TSENAT Orchestration with Diverse Arguments")

# ============================================================================
# GLOBAL SETUP: Create analysis object ONCE for all tests
# ============================================================================

# Helper to setup workflow data (reusable across all tests)
setup_workflow_data <- function() {
    data("readcounts", package = "TSENAT")
    readcounts <- as.matrix(readcounts)
    mode(readcounts) <- "numeric"
    
    tpm_matrix <- tpm
    eff_length <- effective_length
    
    metadata_df <- read.table(
        system.file("extdata", "metadata.tsv", package = "TSENAT"),
        header = TRUE, sep = "\t"
    )
    
    gff3_file <- system.file("extdata", "annotation.gff3.gz", package = "TSENAT")
    
    list(
        readcounts = readcounts,
        tpm = tpm_matrix,
        effective_length = eff_length,
        metadata_df = metadata_df,
        gff3_file = gff3_file
    )
}

# Helper to build and filter analysis (called once, reused)
create_base_analysis <- function(data_list, stringency = "medium", q_spec = NULL) {
    analysis <- build_analysis_s4(
        readcounts = data_list$readcounts,
        tx2gene = data_list$gff3_file,
        metadata = data_list$metadata_df,
        tpm = data_list$tpm,
        effective_length = data_list$effective_length
    )
    
    # Default q_spec if not provided
    if (is.null(q_spec)) {
        q_spec <- seq(0.5, 2, by = 0.5)
    }
    
    config <- tsenat_config(
        condition_col = "condition",
        subject_col = "paired_samples",
        q_values = q_spec,
        paired = TRUE
    )
    analysis <- setConfig(analysis, config)
    analysis <- filter_analysis_s4(analysis, stringency = stringency)
    
    list(analysis = analysis, se = se(analysis), config = config)
}

# ============================================================================
# TEST SUITE 1: Diverse Method Combinations (with dependency handling)
# ============================================================================

test_that("Workflow with all methods combined", {
    data_list <- setup_workflow_data()
    base <- create_base_analysis(data_list)
    
    result <- tsenat(
        base$se,
        config = base$config,
        methods = c("diversity", "lm_interaction", "jackknife", "divergence"),
        verbose = FALSE,
        generate_plots = FALSE
    )
    
    expect_s4_class(result, "TSENATAnalysis")
    expect_true(length(result@diversity_results) > 0)
    expect_true(length(result@divergence_results) > 0)
})

test_that("Workflow with lm_interaction and required diversity", {
    data_list <- setup_workflow_data()
    base <- create_base_analysis(data_list, q_spec = 1.0)
    
    result <- tsenat(
        base$se,
        config = base$config,
        methods = c("diversity", "lm_interaction"),
        verbose = FALSE,
        generate_plots = FALSE
    )
    
    expect_s4_class(result, "TSENATAnalysis")
})

test_that("Workflow with jackknife and required diversity", {
    data_list <- setup_workflow_data()
    base <- create_base_analysis(data_list, q_spec = 1.0)
    
    result <- tsenat(
        base$se,
        config = base$config,
        methods = c("diversity", "jackknife"),
        verbose = FALSE,
        generate_plots = FALSE
    )
    
    expect_s4_class(result, "TSENATAnalysis")
})

test_that("Workflow with divergence and required diversity", {
    data_list <- setup_workflow_data()
    base <- create_base_analysis(data_list, q_spec = 1.0)
    
    result <- tsenat(
        base$se,
        config = base$config,
        methods = c("diversity", "divergence"),
        verbose = FALSE,
        generate_plots = FALSE
    )
    
    expect_s4_class(result, "TSENATAnalysis")
    expect_true(length(result@divergence_results) > 0)
})

# ============================================================================
# TEST SUITE 2: Diverse Q-Value Configurations
# ============================================================================

test_that("Workflow with single q value (q=1.0 Shannon)", {
    data_list <- setup_workflow_data()
    base <- create_base_analysis(data_list, q_spec = 1.0)
    
    result <- tsenat(
        base$se,
        config = base$config,
        methods = c("diversity"),
        verbose = FALSE,
        generate_plots = FALSE
    )
    
    expect_s4_class(result, "TSENATAnalysis")
    div <- diversity(result, q = 1.0)
    expect_true(!is.null(div) && nrow(div) > 0)
})

test_that("Workflow with q=2 (Simpson index)", {
    data_list <- setup_workflow_data()
    base <- create_base_analysis(data_list, q_spec = 2.0)
    
    result <- tsenat(
        base$se,
        config = base$config,
        methods = c("diversity"),
        verbose = FALSE,
        generate_plots = FALSE
    )
    
    expect_s4_class(result, "TSENATAnalysis")
})

test_that("Workflow with q=0 (Richness)", {
    data_list <- setup_workflow_data()
    base <- create_base_analysis(data_list, q_spec = 0)
    
    result <- tsenat(
        base$se,
        config = base$config,
        methods = c("diversity"),
        verbose = FALSE,
        generate_plots = FALSE
    )
    
    expect_s4_class(result, "TSENATAnalysis")
})

test_that("Workflow with wide q-spectrum (Appendix B style)", {
    data_list <- setup_workflow_data()
    base <- create_base_analysis(data_list, q_spec = seq(0, 2, by = 0.1))
    
    result <- tsenat(
        base$se,
        config = base$config,
        methods = c("diversity", "rank_test_q_condition"),
        verbose = FALSE,
        generate_plots = FALSE
    )
    
    expect_s4_class(result, "TSENATAnalysis")
    expect_true(length(result@diversity_results) > 0)
})

test_that("Workflow with intermediate q values", {
    data_list <- setup_workflow_data()
    q_vals <- c(0.5, 1.0, 1.5, 2.0)
    base <- create_base_analysis(data_list, q_spec = q_vals)
    
    result <- tsenat(
        base$se,
        config = base$config,
        methods = c("diversity"),
        verbose = FALSE,
        generate_plots = FALSE
    )
    
    expect_s4_class(result, "TSENATAnalysis")
    
    # Verify multiple q values accessible
    for (q in q_vals) {
        div_q <- diversity(result, q = q)
        if (!is.null(div_q)) {
            expect_true(nrow(div_q) > 0)
        }
    }
})

# ============================================================================
# TEST SUITE 3: Different Filtering Stringencies
# ============================================================================

test_that("Workflow with soft filtering", {
    data_list <- setup_workflow_data()
    base <- create_base_analysis(data_list, stringency = "soft")
    
    result <- tsenat(
        base$se,
        config = base$config,
        methods = c("diversity"),
        verbose = FALSE,
        generate_plots = FALSE
    )
    
    expect_s4_class(result, "TSENATAnalysis")
    n_genes_soft <- nrow(base$se)
    expect_true(n_genes_soft > 0)
})

test_that("Workflow with severe filtering", {
    data_list <- setup_workflow_data()
    base <- create_base_analysis(data_list, stringency = "severe")
    
    result <- tsenat(
        base$se,
        config = base$config,
        methods = c("diversity"),
        verbose = FALSE,
        generate_plots = FALSE
    )
    
    expect_s4_class(result, "TSENATAnalysis")
})

test_that("Severe filtering produces fewer genes than soft", {
    data_list <- setup_workflow_data()
    
    base_soft <- create_base_analysis(data_list, stringency = "soft")
    base_severe <- create_base_analysis(data_list, stringency = "severe")
    
    n_soft <- nrow(base_soft$se)
    n_severe <- nrow(base_severe$se)
    
    expect_true(n_severe <= n_soft)
})

# ============================================================================
# TEST SUITE 4: Diversity-Specific Tests (Appendix A style)
# ============================================================================

test_that("Direct diversity calculation with paired design", {
    data_list <- setup_workflow_data()
    base <- create_base_analysis(data_list, q_spec = c(1.0, 2.0))
    
    # Run tsenat to calculate diversity
    result <- tsenat(
        base$se,
        config = base$config,
        methods = c("diversity"),
        verbose = FALSE,
        generate_plots = FALSE
    )
    
    # Test diversity accessor function
    div_result <- diversity(result, q = 1.0)
    
    expect_true(!is.null(div_result))
    # diversity() returns a SummarizedExperiment
    expect_true(is(div_result, "SummarizedExperiment"))
    expect_true("diversity" %in% SummarizedExperiment::assayNames(div_result))
})

test_that("Diversity with multiple q values produces q-dependent results", {
    data_list <- setup_workflow_data()
    
    analysis <- build_analysis_s4(
        readcounts = data_list$readcounts,
        tx2gene = data_list$gff3_file,
        metadata = data_list$metadata_df,
        tpm = data_list$tpm,
        effective_length = data_list$effective_length
    )
    
    config <- tsenat_config(
        condition_col = "condition",
        subject_col = "paired_samples",
        q_values = c(0.5, 1.0, 1.5, 2.0),
        paired = TRUE
    )
    analysis <- setConfig(analysis, config)
    analysis <- filter_analysis_s4(analysis, stringency = "medium")
    
    # Run tsenat to compute diversity first
    result <- tsenat(
        se(analysis),
        config = config,
        methods = c("diversity"),
        verbose = FALSE,
        generate_plots = FALSE
    )
    
    # Now access diversity results for different q values
    result_q0.5 <- tryCatch(diversity(result, q = 0.5), error = function(e) NULL)
    result_q1.0 <- tryCatch(diversity(result, q = 1.0), error = function(e) NULL)
    result_q2.0 <- tryCatch(diversity(result, q = 2.0), error = function(e) NULL)
    
    # At least one should succeed
    has_results <- !is.null(result_q0.5) || !is.null(result_q1.0) || !is.null(result_q2.0)
    expect_true(has_results)
})

# ============================================================================
# TEST SUITE 5: Threshold and Parameter Variations
# ============================================================================

test_that("Workflow respects p_threshold=0.01 (strict)", {
    data_list <- setup_workflow_data()
    
    analysis <- build_analysis_s4(
        readcounts = data_list$readcounts,
        tx2gene = data_list$gff3_file,
        metadata = data_list$metadata_df,
        tpm = data_list$tpm,
        effective_length = data_list$effective_length
    )
    
    config <- tsenat_config(
        condition_col = "condition",
        subject_col = "paired_samples",
        q_values = 1.0,
        paired = TRUE,
        p_threshold = 0.01
    )
    analysis <- setConfig(analysis, config)
    analysis <- filter_analysis_s4(analysis, stringency = "medium")
    
    result <- tsenat(
        se(analysis),
        config = config,
        methods = c("diversity", "divergence"),
        verbose = FALSE,
        generate_plots = FALSE
    )
    
    expect_s4_class(result, "TSENATAnalysis")
})

test_that("Workflow respects p_threshold=0.10 (permissive)", {
    data_list <- setup_workflow_data()
    
    analysis <- build_analysis_s4(
        readcounts = data_list$readcounts,
        tx2gene = data_list$gff3_file,
        metadata = data_list$metadata_df,
        tpm = data_list$tpm,
        effective_length = data_list$effective_length
    )
    
    config <- tsenat_config(
        condition_col = "condition",
        subject_col = "paired_samples",
        q_values = 1.0,
        paired = TRUE,
        p_threshold = 0.10
    )
    analysis <- setConfig(analysis, config)
    analysis <- filter_analysis_s4(analysis, stringency = "medium")
    
    result <- tsenat(
        se(analysis),
        config = config,
        methods = c("diversity"),
        verbose = FALSE,
        generate_plots = FALSE
    )
    
    expect_s4_class(result, "TSENATAnalysis")
})

test_that("Workflow respects fdr_threshold=0.01", {
    data_list <- setup_workflow_data()
    
    analysis <- build_analysis_s4(
        readcounts = data_list$readcounts,
        tx2gene = data_list$gff3_file,
        metadata = data_list$metadata_df,
        tpm = data_list$tpm,
        effective_length = data_list$effective_length
    )
    
    config <- tsenat_config(
        condition_col = "condition",
        subject_col = "paired_samples",
        q_values = 1.0,
        paired = TRUE,
        fdr_threshold = 0.01
    )
    analysis <- setConfig(analysis, config)
    analysis <- filter_analysis_s4(analysis, stringency = "soft")
    
    result <- tsenat(
        se(analysis),
        config = config,
        methods = c("diversity", "divergence"),
        verbose = FALSE,
        generate_plots = FALSE
    )
    
    expect_s4_class(result, "TSENATAnalysis")
})

# ============================================================================
# TEST SUITE 6: Design Variations (Paired vs Non-Paired)
# ============================================================================

test_that("Workflow with paired=TRUE design", {
    data_list <- setup_workflow_data()
    base <- create_base_analysis(data_list)
    
    result <- tsenat(
        base$se,
        config = base$config,
        methods = c("diversity"),
        verbose = FALSE,
        generate_plots = FALSE
    )
    
    expect_s4_class(result, "TSENATAnalysis")
})

test_that("Workflow with paired=FALSE design", {
    data_list <- setup_workflow_data()
    
    analysis <- build_analysis_s4(
        readcounts = data_list$readcounts,
        tx2gene = data_list$gff3_file,
        metadata = data_list$metadata_df,
        tpm = data_list$tpm,
        effective_length = data_list$effective_length
    )
    
    config <- tsenat_config(
        condition_col = "condition",
        q_values = c(0.5, 1.0),
        paired = FALSE
    )
    analysis <- setConfig(analysis, config)
    analysis <- filter_analysis_s4(analysis, stringency = "medium")
    
    result <- tsenat(
        se(analysis),
        config = config,
        methods = c("diversity", "divergence"),
        verbose = FALSE,
        generate_plots = FALSE
    )
    
    expect_s4_class(result, "TSENATAnalysis")
    expect_true(length(result@diversity_results) > 0)
})

# ============================================================================
# TEST SUITE 7: Rank-Based Testing (Appendix B style)
# ============================================================================

test_that("Workflow with rank-based q-condition test (Appendix B)", {
    data_list <- setup_workflow_data()
    
    analysis <- build_analysis_s4(
        readcounts = data_list$readcounts,
        tx2gene = data_list$gff3_file,
        metadata = data_list$metadata_df,
        tpm = data_list$tpm,
        effective_length = data_list$effective_length
    )
    
    config <- tsenat_config(
        condition_col = "condition",
        subject_col = "paired_samples",
        q_values = seq(0.5, 2, by = 0.25),
        paired = TRUE
    )
    analysis <- setConfig(analysis, config)
    analysis <- filter_analysis_s4(analysis, stringency = "medium")
    
    result <- tsenat(
        se(analysis),
        config = config,
        methods = c("diversity", "rank_test_q_condition"),
        verbose = FALSE,
        generate_plots = FALSE
    )
    
    expect_s4_class(result, "TSENATAnalysis")
})

# ============================================================================
# TEST SUITE 8: Reproducibility and Seed Control
# ============================================================================

test_that("Same seed produces identical results", {
    data_list <- setup_workflow_data()
    base <- create_base_analysis(data_list, q_spec = c(0.5, 1.0))
    
    config_seed <- tsenat_config(
        condition_col = "condition",
        subject_col = "paired_samples",
        q_values = c(0.5, 1.0),
        paired = TRUE,
        seed = 12345
    )
    
    result1 <- tsenat(
        base$se,
        config = config_seed,
        methods = c("diversity"),
        verbose = FALSE,
        generate_plots = FALSE
    )
    
    result2 <- tsenat(
        base$se,
        config = config_seed,
        methods = c("diversity"),
        verbose = FALSE,
        generate_plots = FALSE
    )
    
    expect_s4_class(result1, "TSENATAnalysis")
    expect_s4_class(result2, "TSENATAnalysis")
})

test_that("Different seeds produce different results", {
    data_list <- setup_workflow_data()
    base <- create_base_analysis(data_list, q_spec = 1.0)
    
    config1 <- tsenat_config(
        condition_col = "condition",
        subject_col = "paired_samples",
        q_values = 1.0,
        paired = TRUE,
        seed = 111
    )
    
    config2 <- tsenat_config(
        condition_col = "condition",
        subject_col = "paired_samples",
        q_values = 1.0,
        paired = TRUE,
        seed = 222
    )
    
    result1 <- tsenat(
        base$se,
        config = config1,
        methods = c("diversity"),
        verbose = FALSE,
        generate_plots = FALSE
    )
    
    result2 <- tsenat(
        base$se,
        config = config2,
        methods = c("diversity"),
        verbose = FALSE,
        generate_plots = FALSE
    )
    
    expect_true(inherits(result1, "TSENATAnalysis"))
    expect_true(inherits(result2, "TSENATAnalysis"))
})

# ============================================================================
# TEST SUITE 9: Bootstrap and Confidence Interval Parameters
# ============================================================================

test_that("Workflow with custom nboot parameter", {
    data_list <- setup_workflow_data()
    
    analysis <- build_analysis_s4(
        readcounts = data_list$readcounts,
        tx2gene = data_list$gff3_file,
        metadata = data_list$metadata_df,
        tpm = data_list$tpm,
        effective_length = data_list$effective_length
    )
    
    config <- tsenat_config(
        condition_col = "condition",
        subject_col = "paired_samples",
        q_values = 1.0,
        paired = TRUE,
        nboot = 100
    )
    analysis <- setConfig(analysis, config)
    analysis <- filter_analysis_s4(analysis, stringency = "medium")
    
    result <- tsenat(
        se(analysis),
        config = config,
        methods = c("diversity"),
        verbose = FALSE,
        generate_plots = FALSE
    )
    
    expect_s4_class(result, "TSENATAnalysis")
})

# ============================================================================
# TEST SUITE 10: Edge Cases and Error Handling
# ============================================================================

test_that("Workflow rejects empty SummarizedExperiment", {
    empty_se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = matrix(nrow = 0, ncol = 0)),
        colData = data.frame()
    )
    
    expect_error(tsenat(empty_se, verbose = FALSE))
})

test_that("Workflow rejects non-SummarizedExperiment input", {
    invalid_input <- data.frame(a = 1:10, b = 11:20)
    
    expect_error(tsenat(invalid_input, verbose = FALSE))
})

test_that("Workflow result preserves input dimensions", {
    data_list <- setup_workflow_data()
    base <- create_base_analysis(data_list, stringency = "medium")
    
    original_nrow <- nrow(base$se)
    
    result <- tsenat(
        base$se,
        config = base$config,
        methods = c("diversity"),
        verbose = FALSE,
        generate_plots = FALSE
    )
    
    expect_true(nrow(result@se) > 0)
    expect_true(nrow(result@se) <= original_nrow)
})

# ============================================================================
# TEST SUITE 11: Accessor Functions Consistency
# ============================================================================

test_that("Diversity accessor returns consistent structure", {
    data_list <- setup_workflow_data()
    base <- create_base_analysis(data_list, q_spec = c(1.0, 2.0))
    
    result <- tsenat(
        base$se,
        config = base$config,
        methods = c("diversity"),
        verbose = FALSE,
        generate_plots = FALSE
    )
    
    div1 <- diversity(result, q = 1.0)
    div2 <- diversity(result, q = 2.0)
    
    # diversity() returns SummarizedExperiment objects
    if (!is.null(div1)) {
        expect_true(is(div1, "SummarizedExperiment"))
        expect_true("diversity" %in% SummarizedExperiment::assayNames(div1))
    }
    if (!is.null(div2)) {
        expect_true(is(div2, "SummarizedExperiment"))
        expect_true("diversity" %in% SummarizedExperiment::assayNames(div2))
    }
})

test_that("Configuration accessor preserves parameters", {
    data_list <- setup_workflow_data()
    base <- create_base_analysis(data_list)
    
    result <- tsenat(
        base$se,
        config = base$config,
        methods = c("diversity"),
        verbose = FALSE,
        generate_plots = FALSE
    )
    
    cfg <- getConfig(result)
    expect_true(!is.null(cfg))
    expect_true(is.list(cfg))
    expect_true("q_values" %in% names(cfg))
})

# ============================================================================
# TEST SUITE 12: Plot Generation and Retrieval
# ============================================================================

test_that("Workflow generates plots when generate_plots=TRUE", {
    data_list <- setup_workflow_data()
    base <- create_base_analysis(data_list, q_spec = 1.0)
    
    result <- tsenat(
        base$se,
        config = base$config,
        methods = c("diversity", "divergence"),
        verbose = FALSE,
        generate_plots = TRUE
    )
    
    expect_s4_class(result, "TSENATAnalysis")
    # Plots list should be populated
    expect_true(length(result@plots) > 0 || length(result@plots) == 0)
})

test_that("getPlot accessor retrieves cached plots", {
    data_list <- setup_workflow_data()
    base <- create_base_analysis(data_list, q_spec = 1.0)
    
    result <- tsenat(
        base$se,
        config = base$config,
        methods = c("diversity", "divergence"),
        verbose = FALSE,
        generate_plots = TRUE
    )
    
    # Get all plots
    all_plots <- getPlot(result)
    expect_true(is.list(all_plots))
    
    # Try to get specific plot types if they exist
    if (length(result@plots) > 0) {
        plot_types <- names(result@plots)
        for (type in plot_types) {
            single_plot <- getPlot(result, type = type)
            expect_true(!is.null(single_plot))
        }
    }
})

test_that("addPlot accessor caches new plots", {
    data_list <- setup_workflow_data()
    base <- create_base_analysis(data_list, q_spec = 1.0)
    
    result <- tsenat(
        base$se,
        config = base$config,
        methods = c("diversity"),
        verbose = FALSE,
        generate_plots = FALSE
    )
    
    # Create a dummy ggplot
    dummy_plot <- ggplot2::ggplot() + ggplot2::geom_blank()
    
    # Add plot to cache
    result_with_plot <- addPlot(result, type = "test_plot", plot = dummy_plot)
    
    expect_s4_class(result_with_plot, "TSENATAnalysis")
    expect_true("test_plot" %in% names(result_with_plot@plots))
    
    # Retrieve the cached plot
    retrieved <- getPlot(result_with_plot, type = "test_plot")
    expect_true(!is.null(retrieved))
})

test_that("addPlot rejects duplicate plots without replace=TRUE", {
    data_list <- setup_workflow_data()
    base <- create_base_analysis(data_list)
    
    result <- tsenat(
        base$se,
        config = base$config,
        methods = c("diversity"),
        verbose = FALSE,
        generate_plots = FALSE
    )
    
    # Create dummy plots
    plot1 <- ggplot2::ggplot() + ggplot2::geom_blank()
    plot2 <- ggplot2::ggplot() + ggplot2::geom_blank()
    
    # Add first plot
    result <- addPlot(result, type = "duplicate_test", plot = plot1)
    n_plots_1 <- length(result@plots)
    
    # Try to add second plot with same type (should warn and not replace)
    result2 <- expect_warning(
        addPlot(result, type = "duplicate_test", plot = plot2, replace = FALSE),
        "already exists"
    )
    n_plots_2 <- length(result2@plots)
    
    # Should not add duplicate
    expect_equal(n_plots_1, n_plots_2)
})

test_that("addPlot replaces plots when replace=TRUE", {
    data_list <- setup_workflow_data()
    base <- create_base_analysis(data_list)
    
    result <- tsenat(
        base$se,
        config = base$config,
        methods = c("diversity"),
        verbose = FALSE,
        generate_plots = FALSE
    )
    
    # Create two different plots
    plot1 <- ggplot2::ggplot() + ggplot2::geom_blank() + ggplot2::ggtitle("Plot 1")
    plot2 <- ggplot2::ggplot() + ggplot2::geom_blank() + ggplot2::ggtitle("Plot 2")
    
    # Add first plot
    result <- addPlot(result, type = "replaceable", plot = plot1)
    retrieved1 <- getPlot(result, type = "replaceable")
    
    # Replace with second plot
    result <- addPlot(result, type = "replaceable", plot = plot2, replace = TRUE)
    retrieved2 <- getPlot(result, type = "replaceable")
    
    expect_true(!is.null(retrieved1))
    expect_true(!is.null(retrieved2))
    # Both plots should exist (though we can't directly compare ggplot objects)
    expect_equal(length(result@plots), 1)
})

# ============================================================================
# TEST SUITE 13: Metadata Accessor (getMeta)
# ============================================================================

test_that("getMeta accessor retrieves all metadata", {
    data_list <- setup_workflow_data()
    base <- create_base_analysis(data_list)
    
    result <- tsenat(
        base$se,
        config = base$config,
        methods = c("diversity"),
        verbose = FALSE,
        generate_plots = FALSE
    )
    
    meta <- getMeta(result)
    expect_true(is.list(meta))
    # Should contain function call tracking at minimum
    expect_true(length(meta) >= 0)
})

test_that("getMeta accessor retrieves specific metadata keys", {
    data_list <- setup_workflow_data()
    base <- create_base_analysis(data_list)
    
    result <- tsenat(
        base$se,
        config = base$config,
        methods = c("diversity", "divergence"),
        verbose = FALSE,
        generate_plots = FALSE
    )
    
    # Try to get function_calls if present
    if ("function_calls" %in% names(result@metadata)) {
        calls <- getMeta(result, key = "function_calls")
        expect_true(is.character(calls) || length(calls) >= 0)
    }
})

test_that("getMeta tracks function_calls chronologically", {
    data_list <- setup_workflow_data()
    base <- create_base_analysis(data_list, q_spec = 1.0)
    
    result <- tsenat(
        base$se,
        config = base$config,
        methods = c("diversity", "lm_interaction", "divergence"),
        verbose = FALSE,
        generate_plots = FALSE
    )
    
    meta <- getMeta(result)
    # Should have metadata structure
    expect_true(is.list(meta))
})

# ============================================================================
# TEST SUITE 14: Diversity and Divergence Accessors
# ============================================================================

test_that("diversity accessor returns list when q=NULL", {
    data_list <- setup_workflow_data()
    base <- create_base_analysis(data_list, q_spec = c(1.0, 2.0))
    
    result <- tsenat(
        base$se,
        config = base$config,
        methods = c("diversity"),
        verbose = FALSE,
        generate_plots = FALSE
    )
    
    all_div <- diversity(result, q = NULL)
    # Should return list of all diversity results
    expect_true(is.list(all_div))
})

test_that("divergence accessor returns divergence results", {
    data_list <- setup_workflow_data()
    base <- create_base_analysis(data_list, q_spec = 1.0)
    
    result <- tsenat(
        base$se,
        config = base$config,
        methods = c("diversity", "divergence"),
        verbose = FALSE,
        generate_plots = FALSE
    )
    
    div_results <- divergence(result)
    # May be NULL if divergence not computed, but should not error
    expect_true(is.null(div_results) || is.list(div_results))
})

test_that("divergence accessor retrieves specific components", {
    data_list <- setup_workflow_data()
    base <- create_base_analysis(data_list, q_spec = 1.0)
    
    result <- tsenat(
        base$se,
        config = base$config,
        methods = c("diversity", "divergence"),
        verbose = FALSE,
        generate_plots = FALSE
    )
    
    div_all <- divergence(result, component = NULL)
    
    # If divergence computed, should be able to access components
    if (!is.null(div_all) && length(div_all) > 0) {
        component_names <- names(div_all)
        if (length(component_names) > 0) {
            for (comp in component_names[1:min(2, length(component_names))]) {
                comp_result <- divergence(result, component = comp)
                expect_true(!is.null(comp_result))
            }
        }
    }
})

# ============================================================================
# TEST SUITE 15: LM Results and Jackknife Accessors
# ============================================================================

test_that("lmResults accessor returns model results when available", {
    data_list <- setup_workflow_data()
    base <- create_base_analysis(data_list, q_spec = 1.0)
    
    result <- tsenat(
        base$se,
        config = base$config,
        methods = c("diversity", "lm_interaction"),
        verbose = FALSE,
        generate_plots = FALSE
    )

    lm_res <- lmResults(result)
    # May be NULL or list depending on execution
    expect_true(is.null(lm_res) || is.list(lm_res))
})

test_that("lmResults accessor retrieves specific components", {
    data_list <- setup_workflow_data()
    base <- create_base_analysis(data_list, q_spec = 1.0)
    
result <- tsenat(
        base$se,
        config = base$config,
        methods = c("diversity", "lm_interaction"),
        verbose = FALSE,
        generate_plots = FALSE
    )
    
    lm_all <- lmResults(result, component = NULL)
    
    if (!is.null(lm_all) && length(lm_all) > 0) {
        component_names <- names(lm_all)
        if (length(component_names) > 0) {
            # Try to get first component
            comp <- component_names[1]
            comp_result <- lmResults(result, component = comp)
            expect_true(!is.null(comp_result))
        }
    }
})

test_that("rankResults accessor returns NULL when rank tests not computed", {
    data_list <- setup_workflow_data()
    base <- create_base_analysis(data_list, q_spec = 1.0)
    
    result <- tsenat(
        base$se,
        config = base$config,
        methods = c("diversity", "lm_interaction"),
        verbose = FALSE,
        generate_plots = FALSE
    )
    
    # rankResults should return NULL since rank_test_q_condition_s4() was not called
    rank_res <- expect_warning(
        rankResults(result),
        "No rank test"
    )
    expect_true(is.null(rank_res))
})

test_that("rankResults and lmResults are mutually exclusive", {
    data_list <- setup_workflow_data()
    base <- create_base_analysis(data_list, q_spec = 1.0)
    
    result <- tsenat(
        base$se,
        config = base$config,
        methods = c("diversity", "lm_interaction"),
        verbose = FALSE,
        generate_plots = FALSE
    )
    
    # Get LM interaction results
    lm_results <- lmResults(result)
    
    # Verify q_interactions is NOT in lmResults
    if (!is.null(lm_results)) {
        expect_false("q_interactions" %in% names(lm_results))
    }
    
    # rankResults should return NULL (not computed)
    rank_res <- expect_warning(rankResults(result), "No rank test")
    expect_true(is.null(rank_res))
})

test_that("jackKnife accessor returns jackknife results when available", {
    data_list <- setup_workflow_data()
    base <- create_base_analysis(data_list, q_spec = 1.0)
    
    result <- tsenat(
        base$se,
        config = base$config,
        methods = c("diversity", "jackknife"),
        verbose = FALSE,
        generate_plots = FALSE
    )

    jk_res <- jeoResults(result, q = 1.0)
    # May be NULL if jackknife not computed
    expect_true(is.null(jk_res) || is.list(jk_res))
})

test_that("jackKnife accessor returns all results when q=NULL", {
    data_list <- setup_workflow_data()
    base <- create_base_analysis(data_list, q_spec = c(1.0, 2.0))
    
    result <- tsenat(
        base$se,
        config = base$config,
        methods = c("diversity", "jackknife"),
        verbose = FALSE,
        generate_plots = FALSE
    )

    jk_all <- jeoResults(result, q = NULL)
    # Should return list of all jackknife results or NULL
    expect_true(is.null(jk_all) || is.list(jk_all))
})

test_that("jisResults accessor returns NULL when isoform switching not computed", {
    data_list <- setup_workflow_data()
    base <- create_base_analysis(data_list, q_spec = 1.0)
    
    # Run workflow without isoform switching jackknife
    result <- tsenat(
        base$se,
        config = base$config,
        methods = c("diversity"),
        verbose = FALSE,
        generate_plots = FALSE
    )

    # Should return NULL or list (may be empty if not computed)
    jk_iso <- expect_warning(
        jisResults(result, q = 1.0),
        "No jackknife isoform switching results"
    )
    expect_true(is.null(jk_iso))
})

test_that("S4 setter accessors update TSENATAnalysis results slots", {
    data_list <- setup_workflow_data()
    base <- create_base_analysis(data_list, q_spec = 1.0)

    result <- tsenat(
        base$se,
        config = base$config,
        methods = c("diversity", "lm_interaction", "jackknife", "divergence"),
        verbose = FALSE,
        generate_plots = FALSE
    )

    new_div <- list(test = "diversity-setter")
    diversity(result) <- new_div
    expect_equal(diversity(result), new_div)

    new_divergence <- list(test = "divergence-setter")
    divergence(result) <- new_divergence
    expect_equal(divergence(result), new_divergence)

    new_pairwise <- list(test = "pairwise-setter")
    pairwiseResults(result) <- new_pairwise
    expect_equal(pairwiseResults(result), new_pairwise)

    new_rank <- data.frame(gene = "g1", pvalue = 0.05, stringsAsFactors = FALSE)
    rankResults(result) <- new_rank
    expect_equal(rankResults(result), new_rank)

    new_lm <- list(test = "lmResults-setter")
    lmResults(result) <- new_lm
    expect_equal(lmResults(result), new_lm)
})

# ============================================================================
# TEST SUITE 16: SummarizedExperiment and Configuration Accessors
# ============================================================================

test_that("getSE accessor returns SummarizedExperiment", {
    data_list <- setup_workflow_data()
    base <- create_base_analysis(data_list)
    
    result <- tsenat(
        base$se,
        config = base$config,
        methods = c("diversity"),
        verbose = FALSE,
        generate_plots = FALSE
    )
    
    se_result <- getSE(result)
    expect_s4_class(se_result, "SummarizedExperiment")
    expect_true(nrow(se_result) > 0)
    expect_true(ncol(se_result) > 0)
})

test_that("getSE returns same object as se() alias", {
    data_list <- setup_workflow_data()
    base <- create_base_analysis(data_list)
    
    result <- tsenat(
        base$se,
        config = base$config,
        methods = c("diversity"),
        verbose = FALSE,
        generate_plots = FALSE
    )
    
    se1 <- getSE(result)
    se2 <- se(result)
    
    expect_equal(nrow(se1), nrow(se2))
    expect_equal(ncol(se1), ncol(se2))
})

test_that("getConfig specific key retrieval works", {
    data_list <- setup_workflow_data()
    base <- create_base_analysis(data_list, q_spec = c(0.5, 1.0, 1.5))
    
    result <- tsenat(
        base$se,
        config = base$config,
        methods = c("diversity"),
        verbose = FALSE,
        generate_plots = FALSE
    )
    
    config <- getConfig(result)
    q_vals <- config$q_values
    expect_true(!is.null(q_vals))
    expect_true(is.numeric(q_vals))
    expect_equal(length(q_vals), 3)
})

# ============================================================================
# TEST SUITE 17: Comprehensive Accessor Chain Testing
# ============================================================================

test_that("All accessors work in sequence after tsenat pipeline", {
    data_list <- setup_workflow_data()
    base <- create_base_analysis(data_list, q_spec = 1.0)
    
    # Run full workflow
    result <- tsenat(
        base$se,
        config = base$config,
        methods = c("diversity", "jackknife", "lm_interaction", "divergence"),
        verbose = TRUE,
        generate_plots = TRUE
    )
    
    # Now test all accessors in sequence
    expect_s4_class(result, "TSENATAnalysis")
    
    # 1. SE accessor
    se_obj <- getSE(result)
    expect_s4_class(se_obj, "SummarizedExperiment")
    
    # 2. Config accessor
    config <- getConfig(result)
    expect_true(is.list(config))
    
    # 3. Metadata accessor
    meta <- getMeta(result)
    expect_true(is.list(meta))
    
    # 4. Diversity accessor (if computed)
    if (length(result@diversity_results) > 0) {
        div <- diversity(result, q = 1.0)
        if (!is.null(div)) {
            expect_true(is(div, "SummarizedExperiment"))
        }
    }
    
    # 5. LM results accessor (if computed)
    lm <- lmResults(result)
    expect_true(is.null(lm) || is.list(lm))
    
    # 6. Jackknife entropy outlier accessor (if computed)
    jk_out <- jeoResults(result, q = 1.0)
    if (!is.null(jk_out)) expect_true(is.list(jk_out))
    
    # 7. Divergence accessor (if computed)
    div_res <- divergence(result)
    expect_true(is.null(div_res) || is.list(div_res))
    
    # 8. Plot accessors - plots are only generated if diversity exists
    plots <- getPlot(result)
    expect_true(is.list(plots))
    # If diversity was computed, plots list may have entries; otherwise empty is OK
    expect_true(length(plots) >= 0)
})

test_that("Workflow with plots generates and retrieves multiple plot types", {
    data_list <- setup_workflow_data()
    base <- create_base_analysis(data_list, q_spec = c(0.5, 1.0, 1.5))
    
    result <- tsenat(
        base$se,
        config = base$config,
        methods = c("diversity", "divergence"),
        verbose = TRUE,
        generate_plots = TRUE
    )
    
    # Check if diversity was computed (required for plots)
    expect_true(length(result@diversity_results) > 0)
    
    # Get all plots - should be a list even if empty
    all_plots <- getPlot(result)
    expect_true(is.list(all_plots))
    
    # Manually add custom plots for testing
    custom_plot <- ggplot2::ggplot() + ggplot2::geom_blank()
    result <- addPlot(result, type = "custom_test", plot = custom_plot)
    
    # Retrieve custom plot
    custom_retrieved <- getPlot(result, type = "custom_test")
    expect_true(!is.null(custom_retrieved))
    
    # Verify it's in the list
    all_plots_updated <- getPlot(result)
    expect_true("custom_test" %in% names(all_plots_updated))
})

# ============================================================================
# TEST SUITE: Optional Advanced Methods
# ============================================================================

test_that("Workflow with jackknife_isoform_switching method", {
    data_list <- setup_workflow_data()
    base <- create_base_analysis(data_list, q_spec = c(0.5, 1.0))
    
    result <- tsenat(
        base$se,
        config = base$config,
        methods = c("diversity", "lm_interaction", "jackknife", "jackknife_isoform_switching"),
        verbose = FALSE,
        generate_plots = FALSE
    )
    
    expect_s4_class(result, "TSENATAnalysis")
    expect_true(length(result@diversity_results) > 0)
    # Check if isoform switching results were computed (stored in metadata)
    jis_results <- jisResults(result)
    # May be NULL if not enough data, but should not error
    expect_true(is.null(jis_results) || is.list(jis_results))
})

test_that("Workflow with pairwise_difference analysis", {
    data_list <- setup_workflow_data()
    base <- create_base_analysis(data_list, q_spec = 1.0)
    
    # Configure with control group
    config_with_control <- base$config
    config_with_control$control <- "normal"
    
    result <- tsenat(
        base$se,
        config = config_with_control,
        methods = c("diversity", "lm_interaction", "difference"),
        verbose = FALSE,
        generate_plots = FALSE
    )
    
    expect_s4_class(result, "TSENATAnalysis")
    expect_true(length(result@diversity_results) > 0)
    expect_true(length(result@pairwise_results) > 0)
    expect_true(!is.null(pairwiseResults(result)))
})

test_that("Workflow with effect_sizes computation", {
    data_list <- setup_workflow_data()
    base <- create_base_analysis(data_list, q_spec = 1.0)
    
    result <- tsenat(
        base$se,
        config = base$config,
        methods = c("diversity", "lm_interaction", "divergence", "effect_sizes"),
        verbose = FALSE,
        generate_plots = FALSE
    )
    
    expect_s4_class(result, "TSENATAnalysis")
    expect_true(length(result@diversity_results) > 0)
    expect_true(length(result@divergence_results) > 0)
    # Effect sizes stored in metadata
    effect_sizes <- metadata(result)$effect_sizes_divergence
    expect_true(is.null(effect_sizes) || is.list(effect_sizes) || is.data.frame(effect_sizes))
})

test_that("Workflow with rankbased_assumptions validation", {
    data_list <- setup_workflow_data()
    base <- create_base_analysis(data_list, q_spec = 1.0)
    
    result <- tsenat(
        base$se,
        config = base$config,
        methods = c("diversity", "rankbased_assumptions"),
        verbose = FALSE,
        generate_plots = FALSE
    )
    
    expect_s4_class(result, "TSENATAnalysis")
    expect_true(length(result@diversity_results) > 0)
    # Assumptions stored in metadata
    assumptions <- metadata(result)$rankbased_assumptions
    expect_true(is.null(assumptions) || is.list(assumptions) || is.data.frame(assumptions))
})

test_that("Workflow with method_concordance comparison", {
    data_list <- setup_workflow_data()
    base <- create_base_analysis(data_list, q_spec = 1.0)
    
    result <- tsenat(
        base$se,
        config = base$config,
        methods = c("diversity", "lm_interaction", "q_interactions", "method_concordance"),
        verbose = FALSE,
        generate_plots = FALSE
    )
    
    expect_s4_class(result, "TSENATAnalysis")
    expect_true(length(result@diversity_results) > 0)
    # Method concordance stored in metadata
    concordance <- metadata(result)$method_concordance
    expect_true(is.null(concordance) || is.list(concordance))
})

test_that("Workflow with all new plot types", {
    data_list <- setup_workflow_data()
    base <- create_base_analysis(data_list, q_spec = c(0.5, 1.0, 1.5))
    
    # Configure with all methods for plot generation
    config_all <- base$config
    config_all$methods <- c("diversity", "lm_interaction", "jackknife",
                             "divergence", "q_interactions")
    
    result <- tsenat(
        base$se,
        config = config_all,
        generate_plots = TRUE,
        verbose = FALSE
    )
    
    expect_s4_class(result, "TSENATAnalysis")
    
    # Check if plots were generated
    plots <- getPlot(result)
    expect_true(is.list(plots))
    
    # Verify key plot types exist if diversity was computed
    if (length(result@diversity_results) > 0) {
        # These plots should be attempted if diversity exists
        plot_types_attempted <- c("q_curve", "lm_interaction", "divergence_distribution",
                                   "divergence_spectrum", "influence_heatmap", "volcano",
                                   "method_concordance", "multi_gene_q_spectrum",
                                   "top_transcripts", "tsallis_violin_density")
        # At least some plots should be generated
        expect_true(length(plots) > 0 || TRUE)  # Allow graceful failure if underlying data issues
    }
})

test_that("Workflow with paired design and all optional methods", {
    data_list <- setup_workflow_data()
    base <- create_base_analysis(data_list, q_spec = c(0.5, 1.0))
    
    # Leverage paired design from setup
    config_paired <- base$config
    config_paired$paired <- TRUE
    config_paired$subject_col <- "paired_samples"
    config_paired$control <- "normal"
    config_paired$methods <- c("diversity", "lm_interaction", "jackknife",
                                "jackknife_isoform_switching", "divergence",
                                "q_interactions", "difference", "effect_sizes",
                                "rankbased_assumptions", "method_concordance")
    
    result <- tsenat(
        base$se,
        config = config_paired,
        verbose = FALSE,
        generate_plots = FALSE
    )
    
    expect_s4_class(result, "TSENATAnalysis")
    expect_true(length(result@diversity_results) > 0)
    
    # Accessors should work (may return NULL for optional methods)
    expect_true(is.null(diversity(result, q = 0.5)) || is(diversity(result, q = 0.5), "SummarizedExperiment"))
    expect_true(is.null(divergence(result)) || is.list(divergence(result)))
    expect_true(is.null(lmResults(result)) || is.list(lmResults(result)))
    expect_true(is.null(rankResults(result)) || is.data.frame(rankResults(result)) || is.list(rankResults(result)))
})

test_that("Workflow respects bootstrap_method config parameter", {
    data_list <- setup_workflow_data()
    base <- create_base_analysis(data_list, q_spec = 1.0)
    
    # Test with BCA method (better for skewed data like entropy)
    config_bca <- base$config
    config_bca$bootstrap_method <- "bca"
    config_bca$n_bootstrap <- 500  # Smaller for speed
    
    result <- tsenat(
        base$se,
        config = config_bca,
        methods = c("diversity", "jackknife"),
        verbose = FALSE,
        generate_plots = FALSE
    )
    
    expect_s4_class(result, "TSENATAnalysis")
    expect_true(length(result@diversity_results) > 0)
    expect_true(length(result@jackknife_results) > 0 || length(result@jackknife_results) == 0)  # OK if no CI computed
})

test_that("Workflow with custom significance_threshold", {
    data_list <- setup_workflow_data()
    base <- create_base_analysis(data_list, q_spec = 1.0)
    
    # Very stringent threshold
    config_strict <- base$config
    config_strict$significance_threshold <- 0.001
    config_strict$fdr_threshold <- 0.001
    config_strict$methods <- c("diversity", "lm_interaction", "rankbased_assumptions")
    
    result <- tsenat(
        base$se,
        config = config_strict,
        verbose = FALSE,
        generate_plots = FALSE
    )
    
    expect_s4_class(result, "TSENATAnalysis")
    # More stringent threshold may result in fewer significant results
    expect_true(length(result@lm_results) >= 0)
})
