#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("Usage: fix_indent_multiple_of_4.R <file1> [file2 ...]")
fix_file <- function(file) {
    text <- readLines(file, warn = FALSE)
    ext <- tools::file_ext(file)
    in_r_chunk <- FALSE
    changed <- FALSE
    changed_lines <- 0L
    for (i in seq_along(text)) {
        line <- text[i]
        process <- TRUE
        if (ext %in% c("Rmd", "rmd", "qmd")) {
            # detect start/end of R code chunk
            if (grepl("^```\\s*\\{?\\s*r", line)) {
                in_r_chunk <- TRUE
                next
            }
            if (in_r_chunk && grepl("^```$", line)) {
                in_r_chunk <- FALSE
                next
            }
            process <- in_r_chunk
        }
        if (!process) next
        # Convert leading tabs to 4 spaces, then ensure leading spaces is multiple of 4
        m <- regexpr("^[ \t]*", line)
        lead <- regmatches(line, m)
        if (nchar(lead) == 0) next
        # replace tabs with 4 spaces
        lead2 <- gsub("\t", strrep(" ", 4), lead, perl = TRUE)
        nsp <- nchar(lead2)
        rem <- nsp %% 4L
        if (rem != 0L) {
            add <- 4L - rem
            newlead <- paste0(lead2, strrep(" ", add))
            newline <- sub("^[ \t]*", newlead, line)
            text[i] <- newline
            changed <- TRUE
            changed_lines <- changed_lines + 1L
        } else if (!identical(lead, lead2)) {
            # tabs replaced but already multiple of 4
            newline <- sub("^[ \t]*", lead2, line)
            text[i] <- newline
            changed <- TRUE
            changed_lines <- changed_lines + 1L
        }
    }
    if (changed) {
        writeLines(text, file)
        abs <- normalizePath(file, winslash = "/", mustWork = FALSE)
        cat("Fixed:", abs, sprintf("(%d lines changed)", changed_lines), "\n")
    } else {
        abs <- normalizePath(file, winslash = "/", mustWork = FALSE)
        cat("No change:", abs, "\n")
    }
}

for (f in args) {
    if (file.exists(f)) fix_file(f) else cat("Missing:", f, "\n")
}
