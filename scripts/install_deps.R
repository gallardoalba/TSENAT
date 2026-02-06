#!/usr/bin/env Rscript
# Install all package dependencies (Imports, Depends, LinkingTo, Suggests)
#
# Usage: Rscript scripts/install_deps.R

tryCatch(
    {
        if (!requireNamespace("remotes", quietly = TRUE)) {
            install.packages("remotes", repos = "https://cloud.r-project.org")
        }
        
        cat("Installing package dependencies...\n")
        remotes::install_deps(dependencies = c("Depends", "Imports", "LinkingTo", "Suggests"))
        cat("✓ Dependencies installed successfully\n")
    },
    error = function(e) {
        cat("✗ ERROR:", conditionMessage(e), "\n")
        quit(status = 1)
    }
)
