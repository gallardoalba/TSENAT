#!/usr/bin/env Rscript
tryCatch(
    {
        pkgload::load_all(".")
        if (!requireNamespace("testthat", quietly = TRUE)) stop("testthat required")
        cat("Running tests...\n")
        testthat::test_dir("tests/testthat", reporter = "summary")
    },
    error = function(e) {
        cat("ERROR:", conditionMessage(e), "\n")
        quit(status = 1)
    }
)
