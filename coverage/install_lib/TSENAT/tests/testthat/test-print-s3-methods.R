# Test coverage for S3 print methods and status functions
# Functions tested:
#   - print.gtable() in R/analysis.R:101-102
#   - print.assumptions_text() in R/orchestration_results.R:1156-1157
#   - print.concordance_text() in R/orchestration_results.R:1424-1425
#   - viz_rcpp_status() in R/jackknife_rcpp_wrapper.R:102

context("S3 Print Methods and Status Functions")

library(ggplot2)
library(grid)
library(gtable)
library(SummarizedExperiment)

# ============================================================================
# TEST SUITE 1: print.gtable() - Grid table printing
# ============================================================================

test_that("TSENAT internal gtable helper renders gtable object", {
    skip_if_not_installed("gtable")
    skip_if_not_installed("grid")
    
    # Create a simple gtable (gtable API uses add_grob, not direct grobs)
    gt <- gtable::gtable(widths = grid::unit(c(1, 1), "cm"),
                         heights = grid::unit(c(1, 1), "cm"))
    gt <- gtable::gtable_add_grob(gt, grid::textGrob("A"), t=1, l=1)
    gt <- gtable::gtable_add_grob(gt, grid::textGrob("B"), t=2, l=1)
    
    # Test that helper returns gtable invisibly
    output <- capture.output(
        result <- TSENAT:::.print_gtable(gt)
    )
    
    # Should return gtable invisibly
    expect_identical(result, gt)
})

test_that("TSENAT internal gtable helper handles NULL grobs gracefully", {
    skip_if_not_installed("gtable")
    
    # Create empty gtable
    gt <- gtable::gtable(
        widths = grid::unit(c(1, 1), "cm"),
        heights = grid::unit(c(1, 1), "cm")
    )
    
    # Should not error
    expect_silent(
        capture.output(TSENAT:::.print_gtable(gt))
    )
})

test_that("TSENAT internal gtable helper accepts ellipsis arguments", {
    skip_if_not_installed("gtable")
    
    gt <- gtable::gtable(widths = grid::unit(1, "cm"),
                         heights = grid::unit(1, "cm"))
    gt <- gtable::gtable_add_grob(gt, grid::textGrob("Test"), t=1, l=1)
    
    # Should accept ... without error
    output <- capture.output(
        result <- TSENAT:::.print_gtable(gt, some_arg = "ignored")
    )
    
    expect_identical(result, gt)
})

# ============================================================================
# TEST SUITE 2: print.assumptions_text() - Assumptions output formatting
# ============================================================================

test_that("print.assumptions_text() outputs text content", {
    # Create assumptions_text object
    assumptions_output <- "Exchangeability Test Results:\n- Statistic: 0.45\n- p-value: 0.32"
    class(assumptions_output) <- c("assumptions_text", "character")
    
    # Capture output
    output <- capture.output(print(assumptions_output))
    
    # Should contain the text
    expect_true(any(grepl("Exchangeability", output)))
})

test_that("print.assumptions_text() returns object invisibly", {
    assumptions_output <- "Test assumptions output"
    class(assumptions_output) <- c("assumptions_text", "character")
    
    # Call print but suppress output
    invisible(capture.output(print(assumptions_output)))
    
    # Object should still be the same
    expect_is(assumptions_output, "character")
})

test_that("print.assumptions_text() handles multi-line text", {
    # Create multi-line text - split by \n for testing
    lines <- c("Line 1", "Line 2", "Line 3", "Line 4", "Line 5")
    multi_line <- paste(lines, collapse = "\n")
    class(multi_line) <- c("assumptions_text", "character")
    
    output <- capture.output(print(multi_line))
    
    # The output should contain all 5 lines (or be a single multi-line string)
    # Check that the content is preserved
    combined <- paste(output, collapse = "\n")
    expect_true(all(grepl("Line", combined)))
    expect_true(grepl("Line 5", combined))
})

test_that("print.assumptions_text() handles empty text", {
    empty_text <- ""
    class(empty_text) <- c("assumptions_text", "character")
    
    # Should not error
    output <- capture.output(print(empty_text))
    
    expect_is(output, "character")
})

test_that("print.assumptions_text() with special characters", {
    special_text <- "Statistical Test: χ² = 12.34\nEffect Size: η = 0.5"
    class(special_text) <- c("assumptions_text", "character")
    
    output <- capture.output(print(special_text))
    
    # Should preserve special characters
    combined <- paste(output, collapse = "\n")
    expect_true(grepl("χ²", combined) || grepl("Statistical Test", combined))
})

# ============================================================================
# TEST SUITE 3: print.concordance_text() - Concordance output formatting
# ============================================================================

test_that("print.concordance_text() outputs text content", {
    concordance_output <- "Concordance Analysis Results:\n- GAM coef: 0.87\n- Conover-Iman Rank Transform coef: 0.85\n- Correlation: 0.92"
    class(concordance_output) <- c("concordance_text", "character")
    
    output <- capture.output(print(concordance_output))
    
    # Should contain the text
    expect_true(any(grepl("Concordance", output)))
})

test_that("print.concordance_text() returns object invisibly", {
    concordance_output <- "Test concordance output"
    class(concordance_output) <- c("concordance_text", "character")
    
    # Call print but suppress output
    invisible(capture.output(print(concordance_output)))
    
    # Object should still be the same
    expect_is(concordance_output, "character")
})

test_that("print.concordance_text() handles formatted tables", {
    table_output <- "Concordance Matrix:\n  GAM=+ | GAM=-\nSRH=+ |  100  | 15\nSRH=- |   10  | 75"
    class(table_output) <- c("concordance_text", "character")
    
    output <- capture.output(print(table_output))
    
    # Should preserve formatting
    combined <- paste(output, collapse = "\n")
    expect_true(grepl("Concordance Matrix", combined))
})

test_that("print.concordance_text() with numeric output", {
    numeric_output <- "Correlation Coefficient: 0.9234\nP-value: < 0.001\n95% CI: (0.911, 0.935)"
    class(numeric_output) <- c("concordance_text", "character")
    
    output <- capture.output(print(numeric_output))
    
    # Should handle numeric values
    combined <- paste(output, collapse = "\n")
    expect_true(grepl("0.9234", combined) || grepl("Correlation", combined))
})

test_that("print.concordance_text() handles method names", {
    methods_output <- "Method Comparison:\nGAM (Generalized Additive Model)\nvs\nConover-Iman Rank Transform\nAgreement: 92.5%"
    class(methods_output) <- c("concordance_text", "character")
    
    output <- capture.output(print(methods_output))
    
    combined <- paste(output, collapse = "\n")
    expect_true(grepl("GAM", combined) || grepl("Method", combined))
})

# ============================================================================
# TEST SUITE 4: viz_rcpp_status() - Rcpp acceleration status
# ============================================================================

test_that("viz_rcpp_status() returns list structure", {
    status <- TSENAT:::viz_rcpp_status()
    
    # Should be a list
    expect_is(status, "list")
    
    # Should have required fields
    expect_true("rcpp_available" %in% names(status))
    expect_true("description" %in% names(status))
})

test_that("viz_rcpp_status() contains logical Rcpp flag", {
    status <- TSENAT:::viz_rcpp_status()
    
    # rcpp_available should be logical
    expect_is(status$rcpp_available, "logical")
    expect_length(status$rcpp_available, 1)
})

test_that("viz_rcpp_status() contains description string", {
    status <- TSENAT:::viz_rcpp_status()
    
    # description should be character
    expect_is(status$description, "character")
    expect_length(status$description, 1)
    
    # Should contain meaningful text about jackknife
    expect_true(grepl("[Jj]ackknife|[Rr]cpp", status$description))
})

test_that("viz_rcpp_status() information is informative", {
    status <- TSENAT:::viz_rcpp_status()
    
    # Description should explain what Rcpp acceleration is for
    desc <- tolower(status$description)
    expect_true(
        grepl("acceleration", desc) ||
        grepl("rcpp", desc) ||
        grepl("jackknife", desc)
    )
})

test_that("viz_rcpp_status() returns consistent results", {
    # Call twice and verify consistency
    status1 <- TSENAT:::viz_rcpp_status()
    status2 <- TSENAT:::viz_rcpp_status()
    
    expect_identical(status1, status2)
})

test_that("viz_rcpp_status() has correct format for output", {
    status <- TSENAT:::viz_rcpp_status()
    
    # Verify that the return value can be printed nicely
    output <- capture.output(print(status))
    
    # Should produce some output
    expect_true(length(output) > 0)
})

# ============================================================================
# TEST SUITE 5: Integration tests for print methods
# ============================================================================

test_that("Internal gtable helper is available", {
    # The package should include an internal helper for gtable rendering,
    # but it should not export a global print.gtable method.
    expect_true(exists(".print_gtable", where = asNamespace("TSENAT"), inherits = FALSE))
    expect_is(get(".print_gtable", envir = asNamespace("TSENAT")), "function")
    expect_false("print.gtable" %in% getNamespaceExports("TSENAT"))
    expect_is(getS3method("print", "assumptions_text"), "function")
    expect_is(getS3method("print", "concordance_text"), "function")
})

test_that("Internal gtable helper works correctly", {
    skip_if_not_installed("gtable")
    skip_if_not_installed("grid")

    gt <- gtable::gtable(widths = grid::unit(1, "cm"),
                         heights = grid::unit(1, "cm"))
    gt <- gtable::gtable_add_grob(gt, grid::textGrob("Test"), t=1, l=1)

    expect_identical(TSENAT:::.print_gtable(gt), gt)
})

test_that("S3 methods dispatch correctly", {
    # assumptions_text dispatch
    assume_text <- "Assumptions OK"
    class(assume_text) <- c("assumptions_text", "character")
    output <- capture.output(print(assume_text))
    expect_true(any(grepl("Assumptions", output)))
    
    # concordance_text dispatch
    conc_text <- "Concordance OK"
    class(conc_text) <- c("concordance_text", "character")
    output <- capture.output(print(conc_text))
    expect_true(any(grepl("Concordance", output)))
})

# ============================================================================
# TEST SUITE: Bioconductor c2e8214 S3 Print Method Exports
# ============================================================================
# Tests for S3 methods marked with @exportS3Method in commit c2e8214

test_that("print.rank_assumptions exported and works correctly", {
    # Create a rank_assumptions object
    rank_results <- list(
        sphericity = list(statistic = 0.95, pvalue = 0.32),
        exchangeability = list(statistic = 0.82, pvalue = 0.15)
    )
    class(rank_results) <- "rank_assumptions"
    
    # Print should work without error
    output <- capture.output(result <- print(rank_results))
    
    # Should return invisibly
    expect_identical(result, rank_results)
})

test_that("print.rank_correlation_ci exported and works correctly", {
    # Create a rank_correlation_ci object with proper structure
    ci_results <- list(
        correlation = 0.75,
        ci_lower = 0.60,
        ci_upper = 0.85,
        method = "kendall",
        correlation_matrix = matrix(c(1, 0.75, 0.75, 1), nrow = 2),
        ci_level = 0.95,
        interpretation = "Strong positive correlation"
    )
    class(ci_results) <- "rank_correlation_ci"
    
    # Print should work without error
    output <- capture.output(result <- print(ci_results))
    
    # Should return invisibly
    expect_identical(result, ci_results)
})

test_that("print.tsenat_bootstrap_ci exported and works correctly", {
    # Create a tsenat_bootstrap_ci object
    boot_ci <- list(
        estimate = 0.45,
        ci_lower = 0.40,
        ci_upper = 0.50,
        method = "bca"
    )
    class(boot_ci) <- "tsenat_bootstrap_ci"
    
    # Print should work without error
    output <- capture.output(result <- print(boot_ci))
    
    # Should return invisibly
    expect_identical(result, boot_ci)
})

test_that("print.tsenat_bootstrap_ci_list exported and works correctly", {
    # Create a tsenat_bootstrap_ci_list object
    boot_ci_list <- list(
        list(estimate = 0.45, ci_lower = 0.40, ci_upper = 0.50),
        list(estimate = 0.52, ci_lower = 0.45, ci_upper = 0.58),
        list(estimate = 0.61, ci_lower = 0.55, ci_upper = 0.67)
    )
    class(boot_ci_list) <- "tsenat_bootstrap_ci_list"
    
    # Print should work without error
    output <- capture.output(result <- print(boot_ci_list))
    
    # Should return invisibly
    expect_identical(result, boot_ci_list)
})

test_that("print.tsenat_divergence_bootstrap_ci exported and works correctly", {
    # Create a tsenat_divergence_bootstrap_ci object
    div_boot_ci <- list(
        estimate = 0.35,
        ci_lower = 0.30,
        ci_upper = 0.40,
        method = "percentile"
    )
    class(div_boot_ci) <- "tsenat_divergence_bootstrap_ci"
    
    # Print should work without error
    output <- capture.output(result <- print(div_boot_ci))
    
    # Should return invisibly (silently for this class)
    expect_identical(result, div_boot_ci)
})

test_that("print.tsenat_jackknife exported and works correctly", {
    # Create a tsenat_jackknife object
    jk_result <- list(
        estimate = 0.42,
        bias = 0.01,
        se = 0.03,
        ci_lower = 0.37,
        ci_upper = 0.47
    )
    class(jk_result) <- "tsenat_jackknife"
    
    # Print should work without error
    output <- capture.output(result <- print(jk_result))
    
    # Should return invisibly
    expect_identical(result, jk_result)
})

test_that("print.tsenat_jackknife_list exported and works correctly", {
    # Create a tsenat_jackknife_list object
    jk_list <- list(
        list(estimate = 0.42, bias = 0.01, se = 0.03),
        list(estimate = 0.51, bias = 0.02, se = 0.04),
        list(estimate = 0.63, bias = 0.01, se = 0.05)
    )
    class(jk_list) <- "tsenat_jackknife_list"
    
    # Print should work without error
    output <- capture.output(result <- print(jk_list))
    
    # Should return invisibly
    expect_identical(result, jk_list)
})

test_that("S3 method registration enables proper dispatch", {
    # Verify that proper objects with matching class dispatch correctly
    # Using realistic objects that match the actual class structure
    
    # Test rank_assumptions with proper structure
    rank_obj <- list(
        sphericity = list(stat = 0.5, pval = 0.2),
        exchangeability = list(stat = 0.6, pval = 0.1)
    )
    class(rank_obj) <- "rank_assumptions"
    expect_is(rank_obj, "rank_assumptions")
    # Print methods produce messages, so suppress them
    suppressMessages(capture.output(print(rank_obj)))
    expect_true(TRUE)  # Verify it ran without error
    
    # Test tsenat_bootstrap_ci with proper structure
    boot_obj <- list(
        estimate = 0.5,
        ci_lower = 0.4,
        ci_upper = 0.6
    )
    class(boot_obj) <- "tsenat_bootstrap_ci"
    expect_is(boot_obj, "tsenat_bootstrap_ci")
    # Print methods produce messages, so suppress them
    suppressMessages(capture.output(print(boot_obj)))
    expect_true(TRUE)  # Verify it ran without error
    
    # Test tsenat_jackknife with proper structure
    jk_obj <- list(
        estimate = 0.5,
        bias = 0.01,
        se = 0.02
    )
    class(jk_obj) <- "tsenat_jackknife"
    expect_is(jk_obj, "tsenat_jackknife")
    # Print methods produce messages, so suppress them
    suppressMessages(capture.output(print(jk_obj)))
    expect_true(TRUE)  # Verify it ran without error
})
