# Setup for TSENAT test suite
# ===========================
# This file is automatically sourced by testthat before running tests.
# It loads utility functions and factory functions used across test files.

# Source the test data factory functions
# Use robust path resolution that works from any working directory
factory_file <- NULL

# Method 1: Try direct file in current directory (when run from testthat directory)
if (file.exists("tests-factory.R")) {
  factory_file <- "tests-factory.R"
}

# Method 1b: Try via sys.source location if running under testthat
if (is.null(factory_file)) {
  # When run by testthat, we're in tests/testthat/ working directory
  testthat_dir <- system.file("testthat", package = "TSENAT")
  if (testthat_dir != "" && file.exists(file.path(testthat_dir, "tests-factory.R"))) {
    factory_file <- file.path(testthat_dir, "tests-factory.R")
  }
}

# Method 2: Try relative to package root
if (is.null(factory_file) && file.exists(file.path("tests", "testthat", "tests-factory.R"))) {
  factory_file <- file.path("tests", "testthat", "tests-factory.R")
}

# Method 3: Try from parent directory (handle various working directory scenarios)
if (is.null(factory_file) && file.exists(file.path("..", "..", "tests", "testthat", "tests-factory.R"))) {
  factory_file <- file.path("..", "..", "tests", "testthat", "tests-factory.R")
}

# Method 4: Fallback - walk up the directory tree to find it
if (is.null(factory_file)) {
  pkg_root <- NULL
  tryCatch({
    # Try to find package root via DESCRIPTION file
    root <- getwd()
    while (root != dirname(root)) {
      if (file.exists(file.path(root, "DESCRIPTION"))) {
        pkg_root <<- root  # Found package root
        break
      }
      root <- dirname(root)
    }
  }, error = function(e) NULL)
  
  if (!is.null(pkg_root) && file.exists(file.path(pkg_root, "tests", "testthat", "tests-factory.R"))) {
    factory_file <- file.path(pkg_root, "tests", "testthat", "tests-factory.R")
  }
}

if (!is.null(factory_file) && file.exists(factory_file)) {
  source(factory_file, local = FALSE)
} else {
  warning("Could not find tests-factory.R at: ", paste(c("tests-factory.R", "tests/testthat/tests-factory.R", "../../tests/testthat/tests-factory.R"), collapse = " | "), 
          ". Test factory functions will not be available. Working directory: ", getwd())
}

# Silence package startup messages during tests
suppressPackageStartupMessages({
  library(testthat)
  library(SummarizedExperiment)
})
