#!/usr/bin/env Rscript
# Regenerate documentation from roxygen2 comments
#
# Usage: Rscript scripts/run_roxygen.R

tryCatch(
    {
        if (!requireNamespace("roxygen2", quietly = TRUE)) {
            install.packages("roxygen2", repos = "https://cloud.r-project.org")
        }
        cat("Regenerating documentation from roxygen2 comments...\n")
        roxygen2::roxygenise(".")
        cat("✓ Documentation regenerated successfully\n")
    },
    error = function(e) {
        cat("✗ ERROR:", conditionMessage(e), "\n")
        quit(status = 1)
    }
)
