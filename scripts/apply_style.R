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

# Check styler compliance for R package directory
cat("Checking code style with 4-space indentation...\n\n")

# Use style_pkg with dry.run to check without modifying files
result <- style_pkg(
    indent_by = 4,
    strict = FALSE,
    dry = "on"
)

# Debug: print result structure
cat("Debug - Result class:", class(result), "\n")
cat("Debug - Result structure:\n")
str(result)
cat("\n")

# Check if any files needed styling
# style_pkg returns TRUE if changes were needed, FALSE otherwise
if (isTRUE(result)) {
    cat("ERROR: Code style violations found!\n")
    cat("Please run 'Rscript scripts/apply_style.R' locally to fix styling.\n")
    quit(status = 1)
} else if (is.list(result) && length(result) > 0) {
    cat("ERROR: Code style violations found!\n")
    cat("Please run 'Rscript scripts/apply_style.R' locally to fix styling.\n")
    quit(status = 1)
} else {
    cat("All files are properly styled!\n")
    quit(status = 0)
}
