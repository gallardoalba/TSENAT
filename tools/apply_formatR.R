#!/usr/bin/env Rscript
# Script to apply formatR formatting to R source files
# Works around issues with formatR's file= parameter by reading files first

library(formatR)

if (!requireNamespace("formatR", quietly = TRUE)) {
    install.packages("formatR")
    library(formatR)
}

cat("Applying formatR to R/ files (80 char limit, 4-space indent)...\n")

r_files <- list.files("R", pattern = "\\.R$", full.names = TRUE)
# Skip files with known formatR compatibility issues
skip_files <- c("R/generate_plots.R")
r_files <- setdiff(r_files, skip_files)
failed_files <- character(0)

for (f in r_files) {
    tryCatch({
        # Read file as text instead of using file= parameter
        code_text <- paste(readLines(f, warn = FALSE), collapse = "\n")
        # Format the code - tidy_source returns a list with $text.tidy element
        tidy_result <- tidy_source(text = code_text, width.cutoff = 80, indent = 4, output = FALSE)
        # Extract the formatted code from the list
        formatted_code <- tidy_result$text.tidy
        # Write back to file
        writeLines(formatted_code, f)
        cat(sprintf("✓ %s\n", f))
    }, error = function(e) {
        cat(sprintf("✗ Failed on %s: %s\n", f, conditionMessage(e)))
        failed_files <<- c(failed_files, f)
    })
}

if (length(failed_files) > 0) {
    cat("\nFailed to format the following files:\n")
    for (f in failed_files) cat(sprintf("  - %s\n", f))
    quit(status = 1)
} else {
    cat("\n✓ All files formatted successfully\n")
}
