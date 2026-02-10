#!/usr/bin/env Rscript
# Run package test suite using testthat
#
# Usage: Rscript tools/run_tests.R

tryCatch(
    {
        if (!requireNamespace("pkgload", quietly = TRUE)) {
            install.packages("pkgload", repos = "https://cloud.r-project.org")
        }
        if (!requireNamespace("testthat", quietly = TRUE)) {
            install.packages("testthat", repos = "https://cloud.r-project.org")
        }
        
        cat("Loading package...\n")
        pkgload::load_all(".")
        
        cat("Running tests...\n")
        testthat::test_dir("tests/testthat", reporter = "summary")
        cat("✓ All tests completed\n")
    },
    error = function(e) {
        cat("✗ ERROR:", conditionMessage(e), "\n")
        quit(status = 1)
    }
)
