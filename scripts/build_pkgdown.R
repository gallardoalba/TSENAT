#!/usr/bin/env Rscript
# Build pkgdown documentation site
#
# Usage: Rscript scripts/build_pkgdown.R

tryCatch(
    {
        if (!requireNamespace("pkgdown", quietly = TRUE)) {
            install.packages("pkgdown", repos = "https://cloud.r-project.org")
        }
        
        cat("Building pkgdown site...\n")
        options(warn = -1)  # Suppress spurious stack imbalance warnings
        pkgdown::build_site()
        cat("✓ pkgdown site built successfully in docs/\n")
        invisible()
    },
    error = function(e) {
        cat("✗ ERROR:", conditionMessage(e), "\n")
        quit(status = 1)
    }
)
