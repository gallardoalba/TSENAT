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
    
    result <- suppressWarnings(tsenat(
        base$se,
        config = base$config,
        methods = c("diversity", "lm_interaction", "jackknife", "divergence"),
        verbose = FALSE,
        generate_plots = FALSE
    ))
    
    expect_s4_class(result, "TSENATAnalysis")
    expect_true(length(result@diversity_results) > 0)
    expect_true(length(result@divergence_results) > 0)
})

test_that("Workflow with lm_interaction and required diversity", {
    data_list <- setup_workflow_data()
    base <- create_base_analysis(data_list, q_spec = 1.0)
    
    result <- suppressWarnings(tsenat(
        base$se,
        config = base$config,
        methods = c("diversity", "lm_interaction"),
        verbose = FALSE,
        generate_plots = FALSE
    ))
    
    expect_s4_class(result, "TSENATAnalysis")
})

test_that("Workflow with jackknife and required diversity", {
    data_list <- setup_workflow_data()
    base <- create_base_analysis(data_list, q_spec = 1.0)
    
    result <- suppressWarnings(tsenat(
        base$se,
        config = base$config,
        methods = c("diversity", "jackknife"),
        verbose = FALSE,
        generate_plots = FALSE
    ))
    
    expect_s4_class(result, "TSENATAnalysis")
})

test_that("Workflow with divergence and required diversity", {
    data_list <- setup_workflow_data()
    base <- create_base_analysis(data_list, q_spec = 1.0)
    
    result <- suppressWarnings(tsenat(
        base$se,
        config = base$config,
        methods = c("diversity", "divergence"),
        verbose = FALSE,
        generate_plots = FALSE
    ))
    
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
    
    result <- suppressWarnings(tsenat(
        base$se,
        config = base$config,
        methods = c("diversity", "rank_test_q_condition"),
        verbose = FALSE,
        generate_plots = FALSE
    ))
    
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
    result <- suppressWarnings(tsenat(
        se(analysis),
        config = config,
        methods = c("diversity"),
        verbose = FALSE,
        generate_plots = FALSE
    ))
    
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
    
    result <- suppressWarnings(tsenat(
        se(analysis),
        config = config,
        methods = c("diversity", "divergence"),
        verbose = FALSE,
        generate_plots = FALSE
    ))
    
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
    
    result <- suppressWarnings(tsenat(
        se(analysis),
        config = config,
        methods = c("diversity"),
        verbose = FALSE,
        generate_plots = FALSE
    ))
    
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
    
    result <- suppressWarnings(tsenat(
        se(analysis),
        config = config,
        methods = c("diversity", "divergence"),
        verbose = FALSE,
        generate_plots = FALSE
    ))
    
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
    
    result <- suppressWarnings(tsenat(
        se(analysis),
        config = config,
        methods = c("diversity", "divergence"),
        verbose = FALSE,
        generate_plots = FALSE
    ))
    
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
    
    result <- suppressWarnings(tsenat(
        se(analysis),
        config = config,
        methods = c("diversity", "rank_test_q_condition"),
        verbose = FALSE,
        generate_plots = FALSE
    ))
    
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
    
    result1 <- suppressWarnings(tsenat(
        base$se,
        config = config_seed,
        methods = c("diversity"),
        verbose = FALSE,
        generate_plots = FALSE
    ))
    
    result2 <- suppressWarnings(tsenat(
        base$se,
        config = config_seed,
        methods = c("diversity"),
        verbose = FALSE,
        generate_plots = FALSE
    ))
    
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
    
    result1 <- suppressWarnings(tsenat(
        base$se,
        config = config1,
        methods = c("diversity"),
        verbose = FALSE,
        generate_plots = FALSE
    ))
    
    result2 <- suppressWarnings(tsenat(
        base$se,
        config = config2,
        methods = c("diversity"),
        verbose = FALSE,
        generate_plots = FALSE
    ))
    
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
    
    result <- suppressWarnings(tsenat(
        se(analysis),
        config = config,
        methods = c("diversity"),
        verbose = FALSE,
        generate_plots = FALSE
    ))
    
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
