#!/usr/bin/env Rscript
# Prepare package for release: run full checks and documentation build
#
# Usage: Rscript tools/prepare_release.R

cat("==========================================\n")
cat("TSENAT Release Preparation\n")
cat("==========================================\n\n")

steps_completed <- 0
steps_total <- 7

# Step 1: Check package integrity
cat("[1/", steps_total, "] Running full package checks...\n", sep = "")
result <- system("Rscript tools/check_package.R", ignore.stdout = FALSE, ignore.stderr = FALSE)
if (result == 0) {
    cat("✓ Package checks passed\n")
    steps_completed <- steps_completed + 1
} else {
    cat("✗ Package checks failed\n")
}

# Step 2: Run linters and style checks
cat("\n[2/", steps_total, "] Running linters and style checks...\n", sep = "")
result <- system("Rscript tools/run_linters.R", ignore.stdout = FALSE, ignore.stderr = FALSE)
if (result == 0) {
    cat("✓ Linters completed\n")
    steps_completed <- steps_completed + 1
} else {
    cat("⚠ Linters reported issues or failed\n")
}

# Step 3: Run coverage (generates coverage report)
cat("\n[3/", steps_total, "] Running coverage...\n", sep = "")
result <- system("Rscript tools/run_coverage.R", ignore.stdout = FALSE, ignore.stderr = FALSE)
if (result == 0) {
    cat("✓ Coverage report generated\n")
    steps_completed <- steps_completed + 1
} else {
    cat("⚠ Coverage generation failed or skipped\n")
}

# Step 4: Build vignettes
cat("\n[4/", steps_total, "] Building vignettes...\n", sep = "")
result <- system("Rscript tools/build_vignettes.R", ignore.stdout = FALSE, ignore.stderr = FALSE)
if (result == 0) {
    cat("✓ Vignettes built\n")
    steps_completed <- steps_completed + 1
} else {
    cat("✗ Vignette build failed\n")
}

# Step 5: Build pkgdown site
cat("\n[5/", steps_total, "] Building documentation site...\n", sep = "")
result <- system("Rscript tools/build_pkgdown.R", ignore.stdout = FALSE, ignore.stderr = FALSE)
if (result == 0) {
    cat("✓ Documentation site built\n")
    steps_completed <- steps_completed + 1
} else {
    cat("✗ Documentation build failed\n")
}

# Step 6: Check for uncommitted changes
cat("\n[6/", steps_total, "] Checking git status...\n", sep = "")
tryCatch(
    {
        status_output <- system("git status --porcelain", intern = TRUE)
        if (length(status_output) == 0) {
            cat("✓ Working directory is clean\n")
            steps_completed <- steps_completed + 1
        } else {
            cat("⚠ Uncommitted changes detected:\n")
            for (line in status_output) {
                cat("    ", line, "\n")
            }
            cat("    (Consider committing before release)\n")
            steps_completed <- steps_completed + 1
        }

        # Check for unpushed commits (requires an upstream branch)
        upstream_ok <- system("git rev-parse --abbrev-ref --symbolic-full-name @{u} > /dev/null 2>&1", ignore.stdout = TRUE, ignore.stderr = TRUE) == 0
        if (upstream_ok) {
            unpushed <- tryCatch(system("git rev-list --count @{u}..HEAD", intern = TRUE), error = function(e) NA)
            if (!is.na(unpushed) && as.integer(unpushed) > 0) {
                cat(sprintf("    ⚠ You have %s unpushed commit(s). Consider pushing before tagging/releasing.\n", unpushed))
            } else {
                cat("    ✓ No unpushed commits\n")
            }
        } else {
            cat("    ⚠ No upstream branch configured; cannot check for unpushed commits\n")
        }
    },
    error = function(e) {
        cat("⚠ Could not check git status\n")
    }
)

# Step 7: Verify DESCRIPTION and NEWS
cat("\n[7/", steps_total, "] Verifying metadata files...\n", sep = "")
tryCatch(
    {
        desc <- read.dcf("DESCRIPTION")
        version <- desc[1, "Version"]
        title <- desc[1, "Title"]
        
        cat("    Package: TSENAT\n")
        cat("    Version: ", version, "\n", sep = "")
        cat("    Title: ", title, "\n", sep = "")
        
        if (file.exists("NEWS.md")) {
            news_lines <- readLines("NEWS.md", n = 5)
            if (grepl(version, news_lines[1])) {
                cat("    ✓ NEWS.md is up-to-date with version\n")
            } else {
                cat("    ⚠ NEWS.md may not be updated for this version\n")
            }
        }
        
        cat("✓ Metadata files verified\n")
        steps_completed <- steps_completed + 1
    },
    error = function(e) {
        cat("✗ Error reading metadata:", conditionMessage(e), "\n")
    }
)

# Summary
cat("\n==========================================\n")
cat("Release Preparation Summary\n")
cat("==========================================\n")
cat("Completed: ", steps_completed, "/", steps_total, " steps\n\n", sep = "")

if (steps_completed == steps_total) {
    cat("✓ Package is ready for release!\n")
    cat("\nNext steps:\n")
    cat("1. Review vignettes and docs in docs/\n")
    cat("2. Create git tag: git tag v<version>\n")
    cat("3. Push to remote: git push origin stable --tags\n")
} else {
    cat("⚠ Some preparation steps failed or were skipped.\n")
    cat("  Please review and resolve issues before releasing.\n")
}

cat("\n")
if (steps_completed == steps_total) {
    quit(status = 0)
} else {
    quit(status = 1)
}
