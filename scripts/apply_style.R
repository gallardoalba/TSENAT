#!/usr/bin/env Rscript
# Script to check code styling of R files using the styler package
# Uses 4-space indentation
# Exits with code 1 if files don't meet style standards

# Install styler if needed
if (!requireNamespace("styler", quietly = TRUE)) {
    cat("Installing styler package...\n")
    install.packages("styler")
}

library(styler)

# Command-line: use '--apply' to write style changes (dry = "off").
args <- commandArgs(trailingOnly = TRUE)
apply_fixes <- any(args %in% c("--apply", "-a"))

# Check styler compliance for R package directory
cat("Checking code style with 4-space indentation...\n\n")

# Use style_pkg: dry run by default, or apply fixes when requested
if (apply_fixes) {
    cat("Applying style fixes (dry = off)...\n")
    result <- style_pkg(indent_by = 4, strict = FALSE, dry = "off")
} else {
    result <- style_pkg(indent_by = 4, strict = FALSE, dry = "on")
}

# (No debug printing)

# Check if any files needed styling
# Handle different return types from styler::style_pkg():
# - older versions may return a single logical TRUE/FALSE
# - with dry = "on" it returns a data.frame with a `changed` column
if (isTRUE(result)) {
    cat("ERROR: Code style violations found!\n")
    cat("Please run 'Rscript scripts/apply_style.R' locally to fix styling.\n")
    quit(status = 1)
} else if (is.data.frame(result)) {
    # result is a data.frame with files and a logical `changed` column
    if (any(result$changed, na.rm = TRUE)) {
        changed_files <- result$file[which(result$changed)]
        cat("The following files would be changed by styler:\n")
        cat(paste0(" - ", changed_files, "\n"), sep = "")
        if (!apply_fixes) {
            cat("ERROR: Code style violations found!\n")
            cat("Run 'Rscript scripts/apply_style.R --apply' to apply fixes.\n")
            quit(status = 1)
        } else {
            cat("Applied style fixes to the above files.\n")
            quit(status = 0)
        }
    } else {
        cat("All files are properly styled!\n")
        quit(status = 0)
    }
} else if (is.list(result) && length(result) > 0) {
    # Fallback: treat non-empty lists as potential violations
    cat("ERROR: Code style violations found!\n")
    cat("Please run 'Rscript scripts/apply_style.R' locally to fix styling.\n")
    quit(status = 1)
} else {
    cat("All files are properly styled!\n")
    quit(status = 0)
}
