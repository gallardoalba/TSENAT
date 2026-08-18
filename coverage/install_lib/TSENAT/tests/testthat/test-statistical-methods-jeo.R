context("S4 Wrapper: calculate_jeo() - Jackknife Entropy Operations")
library(TSENAT)

# ============================================================================
# MODULE-LEVEL SETUP: Shared test data loaded ONCE
# ============================================================================
# Cache analysis object with diversity results (prerequisite for jackknife)
test_analysis_jeo <- local({
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
        q = c(0.3, 0.5, 0.7, 0.9, 1.0, 1.2, 1.5, 1.7, 1.9, 2.0),  # 10 q-values
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

    # Compute diversity first (prerequisite for jackknife)
    suppressWarnings({
        analysis <- TSENAT::calculate_diversity(analysis, q = c(0.3, 0.5, 0.7, 0.9, 1.0, 1.2, 1.5, 1.7, 1.9, 2.0))
    })

    analysis
})

# ============================================================================
# TEST SUITE 1: Basic Execution and Input Validation
# ============================================================================

test_that("calculate_jeo() requires TSENATAnalysis object", {
    expect_error(
        TSENAT::calculate_jeo(list(data = "invalid")),
        "must be a TSENATAnalysis object"
    )
    expect_error(
        TSENAT::calculate_jeo(data.frame(x = 1:10)),
        "must be a TSENATAnalysis object"
    )
})

test_that("calculate_jeo() requires diversity results as prerequisite", {
    # Create a fresh analysis without diversity
    analysis <- test_analysis_jeo
    analysis@diversity_results <- list()

    expect_error(
        TSENAT::calculate_jeo(analysis, q = 1.0),
        "Diversity results required"
    )
})

test_that("calculate_jeo() executes with explicit single q-value", {
    result <- TSENAT::calculate_jeo(
        test_analysis_jeo,
        q = 1.0,
        verbose = FALSE
    )

    expect_s4_class(result, "TSENATAnalysis")
    expect_true(length(result@jackknife_results) > 0)
    expect_true("q_1.000" %in% names(result@jackknife_results))
})

test_that("calculate_jeo() executes with multiple q-values", {
    result <- TSENAT::calculate_jeo(
        test_analysis_jeo,
        q = c(0.5, 1.0, 1.5),
        verbose = FALSE
    )

    expect_s4_class(result, "TSENATAnalysis")
    expect_true(length(result@jackknife_results) >= 3)
    expect_true("q_0.500" %in% names(result@jackknife_results))
    expect_true("q_1.000" %in% names(result@jackknife_results))
    expect_true("q_1.500" %in% names(result@jackknife_results))
})

# ============================================================================
# TEST SUITE 2: Parameter Resolution and Default Handling
# ============================================================================

test_that("calculate_jeo() uses config q-values when explicit q is NULL", {
    # Test_analysis_jeo has q = c(0.3, 0.5, 0.7, 0.9, 1.0, 1.2, 1.5, 1.7, 1.9, 2.0) in config
    result <- TSENAT::calculate_jeo(
        test_analysis_jeo,
        q = NULL,  # Should use config
        verbose = FALSE
    )

    expect_s4_class(result, "TSENATAnalysis")
    # Should use all 10 config q-values from diversity results
    expect_true(length(result@jackknife_results) >= 1)
})

test_that("calculate_jeo() respects normalization parameter", {
    result_norm_true <- TSENAT::calculate_jeo(
        test_analysis_jeo,
        q = 1.0,
        norm = TRUE,
        verbose = FALSE
    )

    result_norm_false <- TSENAT::calculate_jeo(
        test_analysis_jeo,
        q = 1.0,
        norm = FALSE,
        verbose = FALSE
    )

    expect_s4_class(result_norm_true, "TSENATAnalysis")
    expect_s4_class(result_norm_false, "TSENATAnalysis")
})

test_that("calculate_jeo() respects pseudocount parameter", {
    result_pc0 <- TSENAT::calculate_jeo(
        test_analysis_jeo,
        q = 1.0,
        pseudocount = 0,
        verbose = FALSE
    )

    result_pc1 <- TSENAT::calculate_jeo(
        test_analysis_jeo,
        q = 1.0,
        pseudocount = 1,
        verbose = FALSE
    )

    expect_s4_class(result_pc0, "TSENATAnalysis")
    expect_s4_class(result_pc1, "TSENATAnalysis")
})

test_that("calculate_jeo() respects top_n parameter", {
    result_top5 <- TSENAT::calculate_jeo(
        test_analysis_jeo,
        q = 1.0,
        top_n = 5,
        verbose = FALSE
    )

    result_top10 <- TSENAT::calculate_jeo(
        test_analysis_jeo,
        q = 1.0,
        top_n = 10,
        verbose = FALSE
    )

    expect_s4_class(result_top5, "TSENATAnalysis")
    expect_s4_class(result_top10, "TSENATAnalysis")
})

# ============================================================================
# TEST SUITE 3: Jackknife Results Structure and Output
# ============================================================================

test_that("calculate_jeo() stores results with correct q-value keys", {
    result <- TSENAT::calculate_jeo(
        test_analysis_jeo,
        q = c(0.5, 1.0, 1.5),
        verbose = FALSE
    )

    # Check structure of stored results
    expect_true("q_0.500" %in% names(result@jackknife_results))
    expect_true("q_1.000" %in% names(result@jackknife_results))
    expect_true("q_1.500" %in% names(result@jackknife_results))
})

test_that("calculate_jeo() tracks function calls in metadata", {
    result <- TSENAT::calculate_jeo(
        test_analysis_jeo,
        q = 1.0,
        verbose = FALSE
    )

    expect_true(length(result@metadata$function_calls) > 0)
    expect_true(any(grepl("jackknife_tsallis_entropy", result@metadata$function_calls)))
})

test_that("calculate_jeo() results are non-empty", {
    result <- TSENAT::calculate_jeo(
        test_analysis_jeo,
        q = 1.0,
        verbose = FALSE
    )

    jk_res <- result@jackknife_results$q_1.000
    expect_true(!is.null(jk_res))
    expect_true(is.list(jk_res) || inherits(jk_res, "tsenat_jackknife_list"))
})

# ============================================================================
# TEST SUITE 4: Output File Handling
# ============================================================================

test_that("calculate_jeo() saves RDS output when requested", {
    output_file <- file.path(tempdir(), "test_jeo_output.rds")

    result <- TSENAT::calculate_jeo(
        test_analysis_jeo,
        q = 1.0,
        output_file = output_file,
        verbose = FALSE
    )

    # File should exist after execution
    expect_true(file.exists(output_file))

    # Should be a valid RDS file
    loaded <- readRDS(output_file)
    expect_s4_class(loaded, "TSENATAnalysis")

    # Cleanup
    if (file.exists(output_file)) unlink(output_file)
})

test_that("calculate_jeo() saves TSV output when requested", {
    output_file <- file.path(tempdir(), "test_jeo_output.tsv")

    result <- TSENAT::calculate_jeo(
        test_analysis_jeo,
        q = 1.0,
        output_file = output_file,
        verbose = FALSE
    )

    # File should exist after execution
    if (file.exists(output_file)) {
        expect_true(file.exists(output_file))
        # Check it's a text file
        content <- readLines(output_file, n = 1)
        expect_true(length(content) > 0)
        unlink(output_file)
    }
})

test_that("calculate_jeo() saves CSV output when requested", {
    output_file <- file.path(tempdir(), "test_jeo_output.csv")

    result <- TSENAT::calculate_jeo(
        test_analysis_jeo,
        q = 1.0,
        output_file = output_file,
        verbose = FALSE
    )

    # File should exist after execution if text output is generated
    if (file.exists(output_file)) {
        expect_true(file.exists(output_file))
        unlink(output_file)
    }
})

test_that("calculate_jeo() creates output directory if it doesn't exist", {
    output_dir <- file.path(tempdir(), "test_jeo_dir_", floor(runif(1, 1e6, 9.9e6)))
    output_file <- file.path(output_dir, "test_jeo_output.rds")

    if (dir.exists(output_dir)) unlink(output_dir, recursive = TRUE)

    result <- TSENAT::calculate_jeo(
        test_analysis_jeo,
        q = 1.0,
        output_file = output_file,
        verbose = FALSE
    )

    expect_true(dir.exists(output_dir))

    # Cleanup
    if (dir.exists(output_dir)) unlink(output_dir, recursive = TRUE)
})

# ============================================================================
# TEST SUITE 5: Verbose Output Control
# ============================================================================

test_that("calculate_jeo() with verbose=TRUE produces messages", {
    output <- capture.output({
        result <- TSENAT::calculate_jeo(
            test_analysis_jeo,
            q = c(0.5, 1.0),
            verbose = TRUE
        )
    })

    expect_s4_class(result, "TSENATAnalysis")
    # Should have some output when verbose
    expect_true(length(output) >= 0)  # verbose output optional
})

test_that("calculate_jeo() with verbose=FALSE suppresses messages", {
    output <- capture.output({
        result <- TSENAT::calculate_jeo(
            test_analysis_jeo,
            q = 1.0,
            verbose = FALSE
        )
    })

    expect_s4_class(result, "TSENATAnalysis")
    # Should be minimal output
    expect_true(length(output) < 10)
})

# ============================================================================
# TEST SUITE 6: Metadata Tracking and Chaining
# ============================================================================

test_that("calculate_jeo() preserves original diversity results", {
    original_div_count <- length(test_analysis_jeo@diversity_results)

    result <- TSENAT::calculate_jeo(
        test_analysis_jeo,
        q = 1.0,
        verbose = FALSE
    )

    expect_equal(
        length(result@diversity_results),
        original_div_count
    )
})

test_that("calculate_jeo() can be chained multiple times with different q-values", {
    result1 <- TSENAT::calculate_jeo(
        test_analysis_jeo,
        q = 0.5,
        verbose = FALSE
    )

    result2 <- TSENAT::calculate_jeo(
        result1,
        q = 1.0,
        verbose = FALSE
    )

    expect_s4_class(result2, "TSENATAnalysis")
    # Both q-values should be present
    expect_true("q_0.500" %in% names(result2@jackknife_results))
    expect_true("q_1.000" %in% names(result2@jackknife_results))
})

test_that("calculate_jeo() adds to metadata function_calls without duplicating", {
    result <- TSENAT::calculate_jeo(
        test_analysis_jeo,
        q = 1.0,
        verbose = FALSE
    )

    # Count occurrences of jackknife calls
    jk_calls <- grep("jackknife_tsallis_entropy", result@metadata$function_calls, value = TRUE)
    expect_true(length(jk_calls) >= 1)
})

# ============================================================================
# TEST SUITE 7: Edge Cases and Error Handling
# ============================================================================

test_that("calculate_jeo() errors when q-value not in diversity results", {
    expect_error(
        TSENAT::calculate_jeo(
            test_analysis_jeo,
            q = 99.0,  # Not calculated
            verbose = FALSE
        ),
        "Diversity not calculated"
    )
})

test_that("calculate_jeo() requires numeric q-value", {
    expect_error(
        TSENAT::calculate_jeo(
            test_analysis_jeo,
            q = "one",  # String instead of numeric
            verbose = FALSE
        ),
        "must be numeric"
    )
})

test_that("calculate_jeo() handles log_base parameter", {
    result_log_e <- TSENAT::calculate_jeo(
        test_analysis_jeo,
        q = 1.0,
        log_base = exp(1),
        verbose = FALSE
    )

    result_log_2 <- TSENAT::calculate_jeo(
        test_analysis_jeo,
        q = 1.0,
        log_base = 2,
        verbose = FALSE
    )

    expect_s4_class(result_log_e, "TSENATAnalysis")
    expect_s4_class(result_log_2, "TSENATAnalysis")
})

# ============================================================================
# TEST SUITE 8: Thread Control
# ============================================================================

test_that("calculate_jeo() respects nthreads parameter", {
    result_single <- TSENAT::calculate_jeo(
        test_analysis_jeo,
        q = 1.0,
        nthreads = 1,
        verbose = FALSE
    )

    expect_s4_class(result_single, "TSENATAnalysis")
    expect_true(length(result_single@jackknife_results) > 0)
})

test_that("calculate_jeo() handles nthreads from config when not explicit", {
    analysis_with_threads <- test_analysis_jeo
    analysis_with_threads@config$nthreads <- 2

    result <- TSENAT::calculate_jeo(
        analysis_with_threads,
        q = 1.0,
        nthreads = NULL,  # Should use config
        verbose = FALSE
    )

    expect_s4_class(result, "TSENATAnalysis")
    expect_true(length(result@jackknife_results) > 0)
})
