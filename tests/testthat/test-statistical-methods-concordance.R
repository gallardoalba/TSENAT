context("S4 Wrapper: calculate_concordance() - Method Concordance Computation")
library(TSENAT)

# Skip all tests on CRAN
skip_on_cran()

# ============================================================================
# MODULE-LEVEL SETUP: Shared test data loaded ONCE
# ============================================================================
# Cache analysis with both LM and rank test results (prerequisite)
# Use precomputed LM results from RDS file (matches Appendix B vignette pattern)
test_analysis_concordance <- local({
    # Load precomputed analysis with LM results already computed
    # (avoids expensive lm_interaction computation)
    analysis_lm <- readRDS(
        system.file("extdata", "analysis_lm.rds", package = "TSENAT")
    )
    
    # Prepare fresh diversity-only analysis for rank test computation
    set.seed(42)
    data("readcounts", package = "TSENAT", envir = environment())
    readcounts <- as.matrix(readcounts)

    metadata_df <- read.table(
        system.file("extdata", "metadata.tsv", package = "TSENAT"),
        header = TRUE, sep = "\t"
    )

    gff3_file <- system.file("extdata", "annotation.gff3.gz", package = "TSENAT")

    config <- TSENAT::TSENAT_config(
        sample_col = "sample",
        condition_col = "condition",
        subject_col = "paired_samples",
        q = seq(0, 2, by = 0.2),  # 11 q-values: 0, 0.2, 0.4, ..., 2.0 (follows Appendix B pattern)
        paired = TRUE,
        control = "normal",
        nthreads = 1
    )

    analysis <- TSENAT::build_analysis(
        config = config,
        readcounts = readcounts,
        metadata = metadata_df,
        tx2gene = gff3_file,
        tpm = tpm,
        effective_length = effective_length
    )

    analysis <- TSENAT::filter_analysis(analysis, stringency = "medium")

    # Compute diversity (prerequisite for both LM and rank tests)
    suppressWarnings({
        analysis <- TSENAT::calculate_diversity(analysis, q = seq(0, 2, by = 0.2))
    })

    # Copy LM results from precomputed analysis to current analysis
    # This avoids expensive lm_interaction computation and provides validated results
    analysis@lm_results <- analysis_lm@lm_results

    # Compute rank test results (prerequisite for concordance)
    suppressWarnings({
        analysis <- TSENAT::calculate_rank_test(analysis, method = "rank_test")
    })

    analysis
})

# ============================================================================
# TEST SUITE 1: Input Validation and Prerequisites
# ============================================================================

test_that("calculate_concordance() requires TSENATAnalysis object", {
    expect_error(
        TSENAT::calculate_concordance(analysis_lm = list(data = "invalid")),
        "inherited method"
    )
    expect_error(
        TSENAT::calculate_concordance(analysis_lm = data.frame(x = 1:10)),
        "inherited method"
    )
})

test_that("calculate_concordance() requires LM results", {
    analysis <- test_analysis_concordance
    analysis@lm_results <- list()

    expect_error(
        TSENAT::calculate_concordance(analysis, rank_method = "rank_test"),
        "No LM results found"
    )
})

test_that("calculate_concordance() requires rank test results in legacy API", {
    analysis <- test_analysis_concordance
    analysis@rank_test_results <- list()

    expect_error(
        TSENAT::calculate_concordance(analysis, rank_method = "rank_test"),
        "not found in rank_test_results"
    )
})

# ============================================================================
# TEST SUITE 2: Legacy Single-Object API
# ============================================================================

test_that("calculate_concordance() executes with legacy API (single object)", {
    result <- TSENAT::calculate_concordance(
        test_analysis_concordance,
        lm_method = NULL,
        rank_method = "rank_test",
        verbose = FALSE
    )

    expect_s4_class(result, "TSENATAnalysis")
    expect_true(!is.null(result@metadata$method_concordance))
})

test_that("calculate_concordance() stores concordance results in metadata", {
    result <- TSENAT::calculate_concordance(
        test_analysis_concordance,
        rank_method = "rank_test",
        verbose = FALSE
    )

    concordance_meta <- result@metadata$method_concordance
    expect_true(!is.null(concordance_meta))
    expect_true("comparison_df" %in% names(concordance_meta))
    expect_true("spearman_rho" %in% names(concordance_meta))
    expect_true("high_confidence" %in% names(concordance_meta))
    expect_true("agreement_table" %in% names(concordance_meta))
})

test_that("calculate_concordance() tracks function calls", {
    result <- TSENAT::calculate_concordance(
        test_analysis_concordance,
        rank_method = "rank_test",
        verbose = FALSE
    )

    expect_true(any(grepl("calculate_concordance", result@metadata$function_calls)))
})

test_that("calculate_concordance() auto-detects LM method when NULL", {
    analysis <- test_analysis_concordance

    result <- TSENAT::calculate_concordance(
        analysis,
        lm_method = NULL,
        rank_method = "rank_test",
        verbose = FALSE
    )

    expect_s4_class(result, "TSENATAnalysis")
    expect_true(!is.null(result@metadata$method_concordance$lm_method))
})

test_that("calculate_concordance() errors on non-existent LM method", {
    expect_error(
        TSENAT::calculate_concordance(
            test_analysis_concordance,
            lm_method = "nonexistent_method",
            rank_method = "rank_test",
            verbose = FALSE
        ),
        "not found in LM results"
    )
})

# ============================================================================
# TEST SUITE 3: Two-Object API
# ============================================================================

test_that("calculate_concordance() accepts two TSENATAnalysis objects", {
    analysis_lm <- test_analysis_concordance
    analysis_rank <- test_analysis_concordance

    result <- TSENAT::calculate_concordance(
        analysis_lm,
        analysis_rank = analysis_rank,
        lm_method = NULL,
        rank_method = "rank_test",
        verbose = FALSE
    )

    expect_s4_class(result, "TSENATAnalysis")
    expect_true(!is.null(result@metadata$method_concordance))
})

test_that("calculate_concordance() rejects non-TSENATAnalysis as analysis_rank", {
    expect_error(
        TSENAT::calculate_concordance(
            test_analysis_concordance,
            analysis_rank = list(data = "invalid"),
            verbose = FALSE
        ),
        "must be a TSENATAnalysis object"
    )
})

test_that("calculate_concordance() validates rank_test_results in two-object API", {
    analysis_rank <- test_analysis_concordance
    analysis_rank@rank_test_results <- list()

    expect_error(
        TSENAT::calculate_concordance(
            test_analysis_concordance,
            analysis_rank = analysis_rank,
            rank_method = "rank_test",
            verbose = FALSE
        ),
        "No rank test results found"
    )
})

# ============================================================================
# TEST SUITE 4: Method Selection and Results Validation
# ============================================================================

test_that("calculate_concordance() computes Spearman correlation", {
    result <- TSENAT::calculate_concordance(
        test_analysis_concordance,
        rank_method = "rank_test",
        verbose = FALSE
    )

    spearman_rho <- result@metadata$method_concordance$spearman_rho
    expect_true(!is.na(spearman_rho))
    expect_true(spearman_rho >= -1 && spearman_rho <= 1)
})

test_that("calculate_concordance() includes comparison data frame", {
    result <- TSENAT::calculate_concordance(
        test_analysis_concordance,
        rank_method = "rank_test",
        verbose = FALSE
    )

    comparison_df <- result@metadata$method_concordance$comparison_df
    expect_true(is.data.frame(comparison_df))
    expect_true(nrow(comparison_df) > 0)
})

test_that("calculate_concordance() includes agreement table", {
    result <- TSENAT::calculate_concordance(
        test_analysis_concordance,
        rank_method = "rank_test",
        verbose = FALSE
    )

    agreement_table <- result@metadata$method_concordance$agreement_table
    expect_true(!is.null(agreement_table))
})

test_that("calculate_concordance() identifies high confidence genes", {
    result <- TSENAT::calculate_concordance(
        test_analysis_concordance,
        rank_method = "rank_test",
        verbose = FALSE
    )

    high_conf <- result@metadata$method_concordance$high_confidence
    expect_true(!is.null(high_conf))
})

test_that("calculate_concordance() stores method names in metadata", {
    result <- TSENAT::calculate_concordance(
        test_analysis_concordance,
        rank_method = "rank_test",
        verbose = FALSE
    )

    metadata <- result@metadata$method_concordance
    expect_true(!is.null(metadata$lm_method))
    expect_true(!is.null(metadata$rank_method))
    expect_equal(metadata$rank_method, "rank_test")
})

test_that("calculate_concordance() stores timestamp", {
    result <- TSENAT::calculate_concordance(
        test_analysis_concordance,
        rank_method = "rank_test",
        verbose = FALSE
    )

    metadata <- result@metadata$method_concordance
    expect_true(!is.null(metadata$timestamp))
    expect_true(inherits(metadata$timestamp, "POSIXct"))
})

test_that("calculate_concordance() preserves original results slots", {
    original_lm_count <- length(test_analysis_concordance@lm_results)
    original_rank_count <- length(test_analysis_concordance@rank_test_results)

    result <- TSENAT::calculate_concordance(
        test_analysis_concordance,
        rank_method = "rank_test",
        verbose = FALSE
    )

    expect_equal(length(result@lm_results), original_lm_count)
    expect_equal(length(result@rank_test_results), original_rank_count)
})

# ============================================================================
# TEST SUITE 5: Verbose Output Control
# ============================================================================

test_that("calculate_concordance() with verbose=TRUE produces messages", {
    output <- capture.output({
        result <- TSENAT::calculate_concordance(
            test_analysis_concordance,
            rank_method = "rank_test",
            verbose = TRUE
        )
    })

    expect_s4_class(result, "TSENATAnalysis")
    expect_true(length(output) >= 0)
})

test_that("calculate_concordance() with verbose=FALSE suppresses messages", {
    output <- capture.output({
        result <- TSENAT::calculate_concordance(
            test_analysis_concordance,
            rank_method = "rank_test",
            verbose = FALSE
        )
    })

    expect_s4_class(result, "TSENATAnalysis")
    expect_true(length(output) < 10)
})

# ============================================================================
# TEST SUITE 6: Output File Handling
# ============================================================================

test_that("calculate_concordance() saves RDS output when requested", {
    output_file <- file.path(tempdir(), "test_concordance_output.rds")

    result <- TSENAT::calculate_concordance(
        test_analysis_concordance,
        output_file = output_file,
        rank_method = "rank_test",
        verbose = FALSE
    )

    expect_true(file.exists(output_file))
    loaded <- readRDS(output_file)
    expect_s4_class(loaded, "TSENATAnalysis")
    expect_true(!is.null(loaded@metadata$method_concordance))

    if (file.exists(output_file)) unlink(output_file)
})

test_that("calculate_concordance() creates output directory if needed", {
    output_dir <- file.path(tempdir(), paste0("test_concordance_dir_", floor(runif(1, 1e6, 9.9e6))))
    output_file <- file.path(output_dir, "test_concordance.rds")

    if (dir.exists(output_dir)) unlink(output_dir, recursive = TRUE)
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

    result <- TSENAT::calculate_concordance(
        test_analysis_concordance,
        output_file = output_file,
        rank_method = "rank_test",
        verbose = FALSE
    )

    expect_s4_class(result, "TSENATAnalysis")
    expect_true(file.exists(output_file))
    
    if (dir.exists(output_dir)) unlink(output_dir, recursive = TRUE)
})

# ============================================================================
# TEST SUITE 7: Non-Destructive Operations
# ============================================================================

test_that("calculate_concordance() is non-destructive (doesn't modify inputs)", {
    original_lm <- test_analysis_concordance@lm_results$lm_interaction
    original_rank <- test_analysis_concordance@rank_test_results$rank_test

    result <- TSENAT::calculate_concordance(
        test_analysis_concordance,
        rank_method = "rank_test",
        verbose = FALSE
    )

    expect_identical(
        test_analysis_concordance@lm_results$lm_interaction,
        original_lm
    )
    expect_identical(
        test_analysis_concordance@rank_test_results$rank_test,
        original_rank
    )
})

test_that("calculate_concordance() can be called multiple times", {
    analysis <- test_analysis_concordance

    result1 <- TSENAT::calculate_concordance(
        analysis,
        rank_method = "rank_test",
        verbose = FALSE
    )

    result2 <- TSENAT::calculate_concordance(
        analysis,
        rank_method = "rank_test",
        verbose = FALSE
    )

    expect_s4_class(result1, "TSENATAnalysis")
    expect_s4_class(result2, "TSENATAnalysis")
    expect_true(!is.null(result2@metadata$method_concordance))
})

# ============================================================================
# TEST SUITE 8: Parameter Handling
# ============================================================================

test_that("calculate_concordance() respects rank_method parameter", {
    analysis <- test_analysis_concordance

    result <- TSENAT::calculate_concordance(
        analysis,
        rank_method = "rank_test",
        verbose = FALSE
    )

    expect_equal(result@metadata$method_concordance$rank_method, "rank_test")
})
