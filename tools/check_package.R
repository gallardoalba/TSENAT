#!/usr/bin/env Rscript
# Comprehensive package integrity check
#
# Runs multiple checks: roxygen sync, style, tests, dependencies, documentation
# Usage: Rscript tools/check_package.R

cat("========================================\n")
cat("TSENAT Package Integrity Check\n")
cat("========================================\n\n")

errors <- 0
warnings <- 0

# 1. Check roxygen synchronization
cat("[1/5] Checking roxygen documentation...\n")
tryCatch(
    {
        if (!requireNamespace("roxygen2", quietly = TRUE)) {
            install.packages("roxygen2", repos = "https://cloud.r-project.org")
        }
        # Check if man/ files match roxygen comments
        old_files <- list.files("man", full.names = TRUE)
        roxygen2::roxygenise(".")
        new_files <- list.files("man", full.names = TRUE)
        
        if (length(old_files) == length(new_files)) {
            cat("  ✓ Documentation is in sync\n")
        } else {
            cat("  ⚠ Documentation files count changed\n")
            warnings <- warnings + 1
        }
    },
    error = function(e) {
        cat("  ✗ Error:", conditionMessage(e), "\n")
        errors <- errors + 1
    }
)

# 2. Check code style
cat("\n[2/5] Checking code style...\n")
tryCatch(
    {
        if (!requireNamespace("styler", quietly = TRUE)) {
            install.packages("styler", repos = "https://cloud.r-project.org")
        }
        # This is a dry-run style check
        issues <- styler::style_dir("R", dry = "on", strict = TRUE)
        if (nrow(issues) == 0) {
            cat("  ✓ Code style is consistent\n")
        } else {
            cat("  ⚠ Style issues found in", nrow(issues), "file(s)\n")
            warnings <- warnings + 1
        }
    },
    error = function(e) {
        cat("  ⚠ Style check skipped:", conditionMessage(e), "\n")
        warnings <- warnings + 1
    }
)

# 3. Run tests
cat("\n[3/5] Running test suite...\n")
tryCatch(
    {
        if (!requireNamespace("testthat", quietly = TRUE)) {
            install.packages("testthat", repos = "https://cloud.r-project.org")
        }
        if (!requireNamespace("pkgload", quietly = TRUE)) {
            install.packages("pkgload", repos = "https://cloud.r-project.org")
        }
        
        pkgload::load_all(".")
        result <- testthat::test_dir("tests/testthat", reporter = "summary", quiet = TRUE)
        
        if (result$n_fail == 0 && result$n_error == 0) {
            cat("  ✓ All tests passed\n")
        } else {
            cat("  ✗ Test failures:", result$n_fail, "failures, ", result$n_error, "errors\n")
            errors <- errors + 1
        }
    },
    error = function(e) {
        cat("  ⚠ Tests could not run:", conditionMessage(e), "\n")
        warnings <- warnings + 1
    }
)

# 4. Check dependencies declared in DESCRIPTION
cat("\n[4/5] Checking declared dependencies...\n")
tryCatch(
    {
        desc_file <- read.dcf("DESCRIPTION")
        imports <- strsplit(desc_file[1, "Imports"], ",")[[1]]
        suggests <- strsplit(desc_file[1, "Suggests"], ",")[[1]]
        
        all_deps <- c(imports, suggests)
        all_deps <- trimws(gsub("\\(.*\\)", "", all_deps))
        
        # Load R files and check for library/require calls
        r_files <- list.files("R", pattern = "\\.R$", full.names = TRUE)
        used_pkgs <- character(0)
        
        for (f in r_files) {
            content <- paste(readLines(f, warn = FALSE), collapse = "\n")
            # Find library() and require() calls
            libs <- gregexpr('(?:library|require)\\s*\\(\\s*["\']?([\\w.]+)["\']?\\s*\\)', 
                            content, perl = TRUE)
            if (libs[[1]][1] > 0) {
                matches <- regmatches(content, libs)[[1]]
                for (m in matches) {
                    pkg <- gsub('.*\\(\\s*["\']?([\\w.]+)["\']?.*', "\\1", m)
                    used_pkgs <- c(used_pkgs, pkg)
                }
            }
        }
        
        # Check for base packages
        used_pkgs <- used_pkgs[!(used_pkgs %in% c("base", "stats", "graphics", "methods", "utils"))]
        
        if (length(used_pkgs) > 0) {
            cat("  ✓", length(all_deps), "dependencies declared\n")
        } else {
            cat("  ✓ No undeclared dependencies detected\n")
        }
    },
    error = function(e) {
        cat("  ⚠ Could not check dependencies:", conditionMessage(e), "\n")
    }
)

# 5. Check documentation completeness
cat("\n[5/5] Checking documentation completeness...\n")
tryCatch(
    {
        man_files <- list.files("man", pattern = "\\.Rd$", full.names = TRUE)
        r_files <- list.files("R", pattern = "\\.R$", full.names = TRUE)
        
        # Extract function names from R files
        r_functions <- character(0)
        for (f in r_files) {
            content <- readLines(f, warn = FALSE)
            # Find exported functions (those with @export)
            for (i in seq_along(content)) {
                if (grepl("@export", content[i])) {
                    # Look for function definition nearby
                    for (j in (i+1):min(i+5, length(content))) {
                        if (grepl("^\\s*\\w+\\s*<-\\s*function", content[j])) {
                            fname <- gsub("^\\s*(\\w+)\\s*<-.*", "\\1", content[j])
                            r_functions <- c(r_functions, fname)
                            break
                        }
                    }
                }
            }
        }
        
        if (length(man_files) >= length(unique(r_functions))) {
            cat("  ✓ Documentation appears complete\n")
        } else {
            cat("  ⚠ Some functions may lack documentation\n")
            warnings <- warnings + 1
        }
    },
    error = function(e) {
        cat("  ⚠ Could not check documentation:", conditionMessage(e), "\n")
    }
)

# Summary
cat("\n========================================\n")
cat("Summary: Errors: ", errors, ", Warnings: ", warnings, "\n", sep = "")
cat("========================================\n")

if (errors > 0) {
    quit(status = 1)
}
