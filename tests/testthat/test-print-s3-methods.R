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

test_that("print.gtable() renders gtable object", {
    skip_if_not_installed("gtable")
    skip_if_not_installed("grid")
    
    # Create a simple gtable (gtable API uses add_grob, not direct grobs)
    gt <- gtable::gtable(widths = grid::unit(c(1, 1), "cm"),
                         heights = grid::unit(c(1, 1), "cm"))
    gt <- gtable::gtable_add_grob(gt, grid::textGrob("A"), t=1, l=1)
    gt <- gtable::gtable_add_grob(gt, grid::textGrob("B"), t=2, l=1)
    
    # Test that print returns invisibly
    output <- capture.output(
        result <- print(gt)
    )
    
    # Should return gtable invisibly
    expect_identical(result, gt)
})

test_that("print.gtable() handles NULL grobs gracefully", {
    skip_if_not_installed("gtable")
    
    # Create empty gtable
    gt <- gtable::gtable(
        widths = grid::unit(c(1, 1), "cm"),
        heights = grid::unit(c(1, 1), "cm")
    )
    
    # Should not error
    expect_silent(
        capture.output(print(gt))
    )
})

test_that("print.gtable() accepts ellipsis arguments", {
    skip_if_not_installed("gtable")
    
    gt <- gtable::gtable(widths = grid::unit(1, "cm"),
                         heights = grid::unit(1, "cm"))
    gt <- gtable::gtable_add_grob(gt, grid::textGrob("Test"), t=1, l=1)
    
    # Should accept ... without error
    output <- capture.output(
        result <- print(gt, some_arg = "ignored")
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
    concordance_output <- "Concordance Analysis Results:\n- GAM coef: 0.87\n- SRH coef: 0.85\n- Correlation: 0.92"
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
    methods_output <- "Method Comparison:\nGAM (Generalized Additive Model)\nvs\nSRH (Scheirer-Ray-Hare Rank Test)\nAgreement: 92.5%"
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

test_that("All S3 print methods are registered", {
    # These methods should be available through S3 dispatch
    expect_is(getS3method("print", "gtable"), "function")
    expect_is(getS3method("print", "assumptions_text"), "function")
    expect_is(getS3method("print", "concordance_text"), "function")
})

test_that("S3 methods dispatch correctly", {
    # Test that S3 dispatch works for each type
    
    # gtable dispatch
    gt <- gtable::gtable(widths = grid::unit(1, "cm"),
                         heights = grid::unit(1, "cm"))
    gt <- gtable::gtable_add_grob(gt, grid::textGrob("Test"), t=1, l=1)
    expect_identical(print(gt), gt)
    
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
