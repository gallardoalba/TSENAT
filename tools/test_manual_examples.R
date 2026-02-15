#!/usr/bin/env Rscript
# Script to extract and run examples from R manual (.Rd) files
# This helps validate that the documented examples actually work

# Try to load the package - use devtools::load_all() if in development, otherwise library()
if (requireNamespace("devtools", quietly = TRUE)) {
    devtools::load_all()
} else {
    library(TSENAT)
}

cat("Testing manual examples from man/*.Rd files...\n\n")

# Create a temporary directory for test outputs
test_dir <- tempdir()
results <- data.frame(
    file = character(0),
    status = character(0),
    message = character(0),
    stringsAsFactors = FALSE
)

# Find all .Rd files in the man directory
rd_files <- list.files("man", pattern = "\\.Rd$", full.names = TRUE)

if (length(rd_files) == 0) {
    cat("No .Rd files found in man/ directory\n")
    quit(status = 0)
}

cat(sprintf("Found %d .Rd files to test\n\n", length(rd_files)))

# Process each Rd file
for (rd_file in rd_files) {
    file_name <- basename(rd_file)
    tryCatch({
        # Parse the Rd file
        rd <- tools::parse_Rd(rd_file)
        
        # Extract examples section - look for the examples tag
        examples_idx <- NULL
        for (i in seq_along(rd)) {
            if (!is.null(attr(rd[[i]], "Rd_tag")) && 
                attr(rd[[i]], "Rd_tag") == "\\examples") {
                examples_idx <- i
                break
            }
        }
        
        if (is.null(examples_idx)) {
            results <- rbind(results, data.frame(
                file = file_name,
                status = "SKIP",
                message = "No examples section",
                stringsAsFactors = FALSE
            ))
            cat(sprintf("⊘ %s - No examples\n", file_name))
            next
        }
        
        # Extract the example code text
        ex_section <- rd[[examples_idx]]
        examples_text <- paste(as.character(ex_section), collapse = "\n")
        
        # Remove leading/trailing whitespace and Rd markup
        examples_text <- gsub("^\\\\examples\\{", "", examples_text)
        examples_text <- gsub("\\}$", "", examples_text)
        examples_text <- trimws(examples_text)
        
        if (nchar(examples_text) == 0) {
            results <- rbind(results, data.frame(
                file = file_name,
                status = "SKIP",
                message = "Empty examples section",
                stringsAsFactors = FALSE
            ))
            cat(sprintf("⊘ %s - Empty examples\n", file_name))
            next
        }
        
        # Create a test file
        test_file <- file.path(test_dir, paste0(gsub("\\.Rd$", "", file_name), "_examples.R"))
        writeLines(examples_text, test_file)
        
        # Try to source the examples in a controlled environment
        result_env <- new.env()
        tryCatch({
            source(test_file, local = result_env, echo = FALSE)
            results <- rbind(results, data.frame(
                file = file_name,
                status = "PASS",
                message = "Examples executed successfully",
                stringsAsFactors = FALSE
            ))
            cat(sprintf("✓ %s\n", file_name))
        }, error = function(e) {
            results <<- rbind(results, data.frame(
                file = file_name,
                status = "FAIL",
                message = conditionMessage(e),
                stringsAsFactors = FALSE
            ))
            cat(sprintf("✗ %s - Error: %s\n", file_name, conditionMessage(e)))
        }, warning = function(w) {
            results <<- rbind(results, data.frame(
                file = file_name,
                status = "WARN",
                message = conditionMessage(w),
                stringsAsFactors = FALSE
            ))
            cat(sprintf("⚠ %s - Warning: %s\n", file_name, conditionMessage(w)))
        })
    }, error = function(e) {
        results <<- rbind(results, data.frame(
            file = file_name,
            status = "ERROR",
            message = sprintf("Failed to parse: %s", conditionMessage(e)),
            stringsAsFactors = FALSE
        ))
        cat(sprintf("✗ %s - Parse error: %s\n", file_name, conditionMessage(e)))
    })
}

# Generate summary
sep <- paste0("\n", paste(rep("=", 70), collapse = ""), "\n")
cat(sep)
cat("Summary:\n")
cat(paste(rep("=", 70), collapse = ""), "\n")

summary_table <- table(results$status)
cat(sprintf("Total files: %d\n", nrow(results)))
for (status in names(summary_table)) {
    cat(sprintf("  %s: %d\n", status, summary_table[status]))
}

# Report failures
failed <- results[results$status %in% c("FAIL", "ERROR"), ]
if (nrow(failed) > 0) {
    cat("\nFailed examples:\n")
    cat(paste(rep("-", 70), collapse = ""), "\n")
    for (i in seq_len(nrow(failed))) {
        cat(sprintf("\n%s (%s):\n%s\n", failed$file[i], failed$status[i], failed$message[i]))
    }
    cat(paste(rep("=", 70), collapse = ""), "\n")
    quit(status = 1)
} else {
    cat("\n✓ All examples passed!\n")
    cat(paste(rep("=", 70), collapse = ""), "\n")
    quit(status = 0)
}
