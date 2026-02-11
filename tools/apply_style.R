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
    cat("Applying style fixes (styler; dry = off)...\n")
    result <- style_pkg(indent_by = 4, strict = TRUE, dry = "off")
    # After styler runs, also run formatR to enforce line width wrapping to 80 cols
    if (!requireNamespace("formatR", quietly = TRUE)) {
        cat("Installing formatR package to reformat long lines...\n")
        install.packages("formatR")
    }
    library(formatR)
    r_files <- list.files("R", pattern = "\\.R$", full.names = TRUE)
    for (f in r_files) {
        # tidy_source will overwrite file with wrapped lines when necessary
        tryCatch({
            tidy_source(file = f, width.cutoff = 80, indent = 4, output = TRUE)
        }, error = function(e) {
            message(sprintf("formatR failed on %s: %s", f, conditionMessage(e)))
        })
    }
} else {
    result <- style_pkg(indent_by = 4, strict = TRUE, dry = "on")
    # Additionally, detect files with lines longer than 80 chars
    r_files <- list.files("R", pattern = "\\.R$", full.names = TRUE)
    long_line_files <- character(0)
    for (f in r_files) {
        lines <- readLines(f, warn = FALSE)
        if (any(nchar(lines, type = "width") > 80)) long_line_files <- c(long_line_files, f)
    }
    if (length(long_line_files) > 0) {
        cat("The following files contain lines longer than 80 characters (consider running with --apply to reformat):\n")
        cat(paste0(" - ", long_line_files, "\n"), sep = "")
    }
}

# After applying fixes (if requested), re-check and ensure there are no
# lines longer than 80 characters and that indentation is 4 spaces.
if (!apply_fixes) {
    # result may be logical or data.frame; handle both
    violations <- FALSE
    if (isTRUE(result)) {
        violations <- TRUE
    } else if (is.data.frame(result) && any(result$changed, na.rm = TRUE)) {
        violations <- TRUE
    } else if (is.list(result) && length(result) > 0) {
        violations <- TRUE
    }

    if (violations) {
        cat("ERROR: Code style violations found!\n")
        cat("Run 'Rscript tools/apply_style.R --apply' to apply fixes, which will:\n")
        cat(" - reformat code with styler using 4-space indentation\n")
        cat(" - wrap long lines to 80 characters with formatR\n")
        quit(status = 1)
    }

    # Check for long lines
    long_line_files <- character(0)
    r_files <- list.files("R", pattern = "\\.R$", full.names = TRUE)
    for (f in r_files) {
        lines <- readLines(f, warn = FALSE)
        if (any(nchar(lines, type = "width") > 80)) long_line_files <- c(long_line_files, f)
    }
    if (length(long_line_files) > 0) {
        cat("The following files contain lines longer than 80 characters:\n")
        cat(paste0(" - ", long_line_files, "\n"), sep = "")
        cat("Consider running with --apply to reformat lines.\n")
        quit(status = 1)
    }

    cat("All files are properly styled and under 80 characters per line!\n")
    quit(status = 0)
} else {
    # We applied fixes; ensure styler and formatR rewrote files and line lengths are OK
    cat("Applied styler fixes; running formatR to wrap long lines to 80 chars...\n")
    if (!requireNamespace("formatR", quietly = TRUE)) {
        cat("Installing formatR package to reformat long lines...\n")
        install.packages("formatR")
    }
    library(formatR)

    r_files <- list.files("R", pattern = "\\.R$", full.names = TRUE)
    problematic <- character(0)
    for (f in r_files) {
        # tidy_source with output = FALSE returns a list containing text.tidy
        ok <- tryCatch({
            res <- tidy_source(file = f, width.cutoff = 80, indent = 4, output = FALSE)
            if (!is.null(res$text.tidy)) {
                writeLines(res$text.tidy, con = f)
                TRUE
            } else {
                FALSE
            }
        }, error = function(e) {
            message(sprintf("formatR failed on %s: %s", f, conditionMessage(e)))
            FALSE
        })
        if (!ok) problematic <- c(problematic, f)
    }

    # Re-run a quick length check
    long_line_files <- character(0)
    for (f in r_files) {
        lines <- readLines(f, warn = FALSE)
        if (any(nchar(lines, type = "width") > 80)) long_line_files <- c(long_line_files, f)
    }

    if (length(problematic) > 0 || length(long_line_files) > 0) {
        cat("ERROR: Some files could not be fully reformatted or still have long lines:\n")
        if (length(problematic) > 0) {
            cat(paste0(" - formatR failed on: ", paste(problematic, collapse = ", "), "\n"))
        }
        if (length(long_line_files) > 0) {
            cat(paste0(" - still long lines in: ", paste(long_line_files, collapse = ", "), "\n"))
        }
        quit(status = 1)
    }

    cat("All R files reformatted with 4-space indentation and wrapped to 80 characters.\n")
    quit(status = 0)
}
