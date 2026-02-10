#!/usr/bin/env Rscript
# Run BiocCheck validation against Bioconductor standards
#
# Usage: Rscript tools/run_bioccheck.R

tryCatch(
    {
        options(repos = c(CRAN = "https://cloud.r-project.org"))
        
        if (!requireNamespace("BiocManager", quietly = TRUE)) {
            install.packages("BiocManager", repos = "https://cloud.r-project.org")
        }
        if (!requireNamespace("BiocCheck", quietly = TRUE)) {
            BiocManager::install("BiocCheck")
        }
        
        cat("Running BiocCheck...\n")
        library(BiocCheck)
        BiocCheck::BiocCheck(".")
        cat("✓ BiocCheck completed\n")
    },
    error = function(e) {
        cat("✗ ERROR:", conditionMessage(e), "\n")
        quit(status = 1)
    }
)
