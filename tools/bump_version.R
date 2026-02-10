#!/usr/bin/env Rscript
# Bump package version in DESCRIPTION, NEWS.md, and _pkgdown.yml
#
# Usage: Rscript scripts/bump_version.R [major|minor|patch]
#        Default: patch

args <- commandArgs(trailingOnly = TRUE)
bump_type <- if (length(args) > 0) args[1] else "patch"

if (!(bump_type %in% c("major", "minor", "patch"))) {
    cat("✗ Invalid bump type:", bump_type, "\n")
    cat("  Usage: Rscript scripts/bump_version.R [major|minor|patch]\n")
    quit(status = 1)
}

tryCatch(
    {
        # Read current version from DESCRIPTION
        desc <- read.dcf("DESCRIPTION")
        current_version <- desc[1, "Version"]
        
        # Parse version
        version_parts <- as.numeric(strsplit(current_version, "\\.")[[1]])
        
        # Bump version
        if (bump_type == "major") {
            version_parts[1] <- version_parts[1] + 1
            version_parts[2] <- 0
            version_parts[3] <- 0
        } else if (bump_type == "minor") {
            version_parts[2] <- version_parts[2] + 1
            version_parts[3] <- 0
        } else if (bump_type == "patch") {
            version_parts[3] <- version_parts[3] + 1
        }
        
        new_version <- paste(version_parts, collapse = ".")
        
        cat("Bumping version:", current_version, "→", new_version, "\n\n")
        
        # Update DESCRIPTION
        desc[1, "Version"] <- new_version
        write.dcf(desc, file = "DESCRIPTION")
        cat("✓ Updated DESCRIPTION\n")
        
        # Update _pkgdown.yml if it exists
        if (file.exists("_pkgdown.yml")) {
            pkgdown_content <- readLines("_pkgdown.yml")
            pkgdown_content <- gsub(
                paste0("version: ", current_version),
                paste0("version: ", new_version),
                pkgdown_content
            )
            writeLines(pkgdown_content, "_pkgdown.yml")
            cat("✓ Updated _pkgdown.yml\n")
        }
        
        # Add entry to NEWS.md if it exists
        if (file.exists("NEWS.md")) {
            news_content <- readLines("NEWS.md")
            new_entry <- c(
                paste("# TSENAT", new_version),
                paste("*", format(Sys.Date(), "%Y-%m-%d")),
                "",
                "- (update details here)",
                ""
            )
            writeLines(c(new_entry, news_content), "NEWS.md")
            cat("✓ Updated NEWS.md\n")
        }
        
        cat("\n✓ Version bumped to", new_version, "\n")
        cat("  Remember to update NEWS.md with changes and commit!\n")
    },
    error = function(e) {
        cat("✗ ERROR:", conditionMessage(e), "\n")
        quit(status = 1)
    }
)
