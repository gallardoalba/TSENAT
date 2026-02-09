#!/usr/bin/env Rscript
# Check that R source files are styled. ONLY styles files under `R/`.

if (!requireNamespace("styler", quietly = TRUE)) {
  install.packages("styler", repos = "https://cloud.r-project.org")
}
library(styler)

# Find all R files in the package's R/ directory
r_files <- list.files("R", pattern = "\\.R$", full.names = TRUE)
if (length(r_files) == 0) {
  cat("No R files to style\n")
  quit(status = 0)
}

for (f in r_files) {
  tryCatch({
    txt <- readLines(f, warn = FALSE)
    styled <- styler::style_text(paste(txt, collapse = "\n"), style = styler::tidyverse_style(indent_by = 4))
    # style_text returns a character vector; join and split to preserve line endings
    new_txt <- unlist(strsplit(styled, "\n", fixed = TRUE))
    if (!identical(txt, new_txt)) {
      writeLines(new_txt, con = f)
    }
  }, error = function(e) {
    message("styler failed for file: ", f, " — ", conditionMessage(e))
    # Fallback: try formatR::tidy_file to at least format the file
    tryCatch({
      if (requireNamespace("formatR", quietly = TRUE)) {
        formatR::tidy_file(f, indent = 4, arrow = TRUE)
        message("formatR::tidy_file applied to: ", f)
      }
    }, error = function(e2) {
      message("formatR fallback also failed for file: ", f, " — ", conditionMessage(e2))
    })
  })
}

# If git is available, fail if there are unstaged/uncommitted changes introduced by styler
res <- system("git status --porcelain", intern = TRUE)
if (length(res) > 0 && nzchar(paste(res, collapse = "\n"))) {
  cat("Styler reformatted files; please run `Rscript scripts/check_style.R` and commit the changes before merging.\n")
  system("git --no-pager diff")
  quit(status = 1)
}

cat("No style changes\n")
# Check that code is styled. If styling would change files, exit with failure.
if (!requireNamespace("styler", quietly = TRUE)) {
  install.packages("styler", repos = "https://cloud.r-project.org")
}
library(styler)

# Run styler with 4-space indentation to format in-place.
# Use a safe wrapper: try default styling, and on error retry using token-level
# styling which avoids processing roxygen/knitr examples.
style_safe <- function(path) {
  style_call <- function(scope = NULL) {
    if (is.null(scope)) {
      styler::style_dir(path = path, style = styler::tidyverse_style(indent_by = 4))
    } else {
      styler::style_dir(path = path, scope = scope, style = styler::tidyverse_style(indent_by = 4))
    }
  }

  tryCatch(
    style_call(NULL),
    error = function(e) {
      message("styler failed for '", path, "': ", conditionMessage(e))
      message("Retrying with token-level styling to avoid roxygen/knitr processing.")
      tryCatch(
        style_call("tokens"),
        error = function(e2) {
          message("Token-level styler retry failed for '", path, "': ", conditionMessage(e2))
        }
      )
    }
  )
}

# Style R source files, tests and scripts (safe wrapper handles errors)
style_safe("R")
if (dir.exists("tests")) style_safe("tests")
if (dir.exists("scripts")) style_safe("scripts")

# If git is available, fail if there are unstaged/uncommitted changes introduced by styler
res <- system("git status --porcelain", intern = TRUE)
if (length(res) > 0 && nzchar(paste(res, collapse = "\n"))) {
  cat("Styler reformatted files; please run `Rscript scripts/check_style.R` and commit the changes before merging.\n")
  system("git --no-pager diff")
  quit(status = 1)
}
cat("No style changes\n")
# Check that code is styled. If styling would change files, exit with failure.
# Ensure `styler` is available
if (!requireNamespace("styler", quietly = TRUE)) install.packages("styler", repos = "https://cloud.r-project.org")
# Load styler so internal helpers are registered
library(styler)
# Run styler with 4-space indentation to format in-place
# Use a safe wrapper: try default styling, and on error retry using token-level
# styling which avoids processing roxygen/knitr examples.
style_safe <- function(path) {
    style_call <- function(scope = NULL) {
        if (is.null(scope)) {
            styler::style_dir(
                path = path,
                style = function(...) styler::tidyverse_style(..., indent_by = 4)
            )
        } else {
            styler::style_dir(
                path = path,
                scope = scope,
                style = function(...) styler::tidyverse_style(..., indent_by = 4)
            )
        }
    }
    tryCatch(
        style_call(NULL),
        error = function(e) {
            message("styler failed for '", path, "': ", conditionMessage(e))
            message("Retrying with token-level styling to avoid roxygen/knitr processing.")
            tryCatch(
                style_call("tokens"),
                error = function(e2) {
                    message("Token-level styler retry failed for '", path, "': ", conditionMessage(e2))
                }
            )
        }
    )
}

# Style R source files, tests and scripts (safe wrapper handles errors)
style_safe("R")
if (dir.exists("tests")) style_safe("tests")
if (dir.exists("scripts")) style_safe("scripts")
# If git is available, fail if there are unstaged/uncommitted changes introduced
# by styler
res <- system("git status --porcelain", intern = TRUE)
if (length(res) > 0 && nzchar(paste(res, collapse = "\n"))) {
    cat("Styler reformatted files; please run `Rscript scripts/check_style.R` and commit the changes before merging.\n")
    system("git --no-pager diff")
    quit(status = 1)
}
cat("No style changes\n")
