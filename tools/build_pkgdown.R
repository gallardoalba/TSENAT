#!/usr/bin/env Rscript
# Build pkgdown documentation site
#
# Usage: Rscript scripts/build_pkgdown.R

tryCatch(
    {
        if (!requireNamespace("pkgdown", quietly = TRUE)) {
            install.packages("pkgdown", repos = "https://cloud.r-project.org")
        }

        cat("Checking whether pkgdown site needs rebuild...\n")

        # Define source paths to consider for rebuild triggers
        src_paths <- c("R", "man", "vignettes", "inst", "DESCRIPTION", "README.md", "pkgdown.yml", "NEWS.md")
        src_files <- unlist(lapply(src_paths, function(p) {
            if (file.exists(p)) {
                if (dir.exists(p)) list.files(p, recursive = TRUE, full.names = TRUE) else p
            } else character(0)
        }))
        src_files <- src_files[file.exists(src_files)]

        src_mtime <- if (length(src_files) > 0) max(file.info(src_files)$mtime, na.rm = TRUE) else as.POSIXct(0, origin = "1970-01-01")

        docs_files <- if (dir.exists("docs")) list.files("docs", recursive = TRUE, full.names = TRUE) else character(0)
        docs_mtime <- if (length(docs_files) > 0) max(file.info(docs_files)$mtime, na.rm = TRUE) else as.POSIXct(0, origin = "1970-01-01")

        # Use git information when available to detect uncommitted changes or commit mismatch
        git_available <- system("git rev-parse --is-inside-work-tree > /dev/null 2>&1", ignore.stdout = TRUE) == 0

        rebuild <- FALSE
        if (!dir.exists("docs")) {
            cat("docs/ missing -> will build\n")
            rebuild <- TRUE
        } else if (length(src_files) == 0) {
            cat("No source files found; skipping build\n")
            rebuild <- FALSE
        } else if (src_mtime > docs_mtime) {
            cat("Source files modified after docs/ -> will build\n")
            rebuild <- TRUE
        } else if (git_available) {
            git_status <- system("git status --porcelain", intern = TRUE)
            if (length(git_status) > 0) {
                cat("Uncommitted git changes detected -> will build\n")
                rebuild <- TRUE
            } else {
                commit <- system("git rev-parse HEAD", intern = TRUE)
                commit_file <- file.path("docs", ".built_from_commit")
                if (!file.exists(commit_file) || !identical(readLines(commit_file, warn = FALSE), commit)) {
                    cat("Repository commit differs from last build -> will build\n")
                    rebuild <- TRUE
                } else {
                    cat("Docs appear up-to-date with current commit -> skipping build\n")
                }
            }
        } else {
            cat("No git repository detected and docs are up-to-date -> skipping build\n")
        }

        if (!rebuild) {
            cat("No rebuild necessary; exiting.\n")
            quit(status = 0)
        }

        # Remove docs/ folder entirely if it exists
        if (dir.exists("docs")) {
            cat("Removing existing docs/ folder...\n")
            system2("chmod", c("-R", "u+w", "docs"))  # Make files writable
            system2("rm", c("-rf", "docs"))
        }

        cat("Building pkgdown site...\n")
        options(warn = -1)  # Suppress spurious stack imbalance warnings
        pkgdown::build_site()

        # Record current commit so future runs can skip if unchanged
        if (git_available) {
            commit <- system("git rev-parse HEAD", intern = TRUE)
            dir.create("docs", showWarnings = FALSE, recursive = TRUE)
            writeLines(commit, file.path("docs", ".built_from_commit"))
        }

        cat("✓ pkgdown site built successfully in docs/\n")
        invisible()
    },
    error = function(e) {
        cat("✗ ERROR:", conditionMessage(e), "\n")
        quit(status = 1)
    }
)
