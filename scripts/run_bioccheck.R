#!/usr/bin/env Rscript
tryCatch(
    {
        options(repos = c(CRAN = "https://cloud.r-project.org"))
        if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
        if (!requireNamespace("BiocCheck", quietly = TRUE)) BiocManager::install("BiocCheck")
        library(BiocCheck)
        cat("Running BiocCheck...\n")
        BiocCheck::BiocCheck(".")
        cat("BiocCheck finished.\n")
    },
    error = function(e) {
        cat("ERROR:", conditionMessage(e), "\n")
        quit(status = 1)
    }
)
